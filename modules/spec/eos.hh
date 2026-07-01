#ifndef HYDRO_EOS_HH
#define HYDRO_EOS_HH

#include "types.hh"

#include <flecsi/util/serialize.hh>

#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos.hpp>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace eos {

/**
 * maximum number of different materials
 */
constexpr std::size_t kMaxMaterials = 10;

struct eos_wrapper {
  using eos_t = singularity::Variant<singularity::IdealGas,
    singularity::SpinerEOSDependsRhoSie,
    singularity::Gruneisen,
    singularity::JWL>;

  eos_wrapper() = default;

  template<class M>
  static constexpr bool is_supported_eos_v =
    std::same_as<std::remove_cvref_t<M>, singularity::IdealGas> ||
    std::same_as<std::remove_cvref_t<M>, singularity::SpinerEOSDependsRhoSie> ||
    std::same_as<std::remove_cvref_t<M>, singularity::Gruneisen> ||
    std::same_as<std::remove_cvref_t<M>, singularity::JWL>;

  template<typename M>
    requires is_supported_eos_v<M>
  eos_wrapper(M const & m) : model_(m), lambda_(model_.nlambda()) {}

  static eos_wrapper from_host_device(eos_wrapper const & host) {
    eos_wrapper dev = host;
    dev.model_.GetOnDevice();
    return dev;
  }

  void finalize_device() {
    model_.Finalize();
  }

  /*!
    Return pressure as a function of density and specific internal energy.
   */
  FLECSI_INLINE_TARGET double pRhoSie(double r, double e) const {
    return model_.PressureFromDensityInternalEnergy(r, e);
  }

  /*!
    Return specific internal energy as a function of density and temperature.
   */
  FLECSI_INLINE_TARGET double eRhoT(double r, double t) const {
    return model_.InternalEnergyFromDensityTemperature(r, t);
  }

  /*!
    Return temperature as a function of density and internal specific energy.
   */
  FLECSI_INLINE_TARGET double tRhoSie(double r, double e) const {
    return model_.TemperatureFromDensityInternalEnergy(r, e);
  }

  /*!
    Return sound speed as a function of density and specific internal energy.
   */
  FLECSI_INLINE_TARGET double cRhoSie(double r, double e) const {
    const double b = model_.BulkModulusFromDensityInternalEnergy(r, e);
    return std::sqrt(singularity::robust::ratio(b, r));
  }

private:
  friend flecsi::util::serial::traits<eos_wrapper>;
  template<class Exec>
  friend struct eos_dispatch;

  eos_wrapper(eos_t const & eos, std::vector<double> const & l)
    : model_(eos), lambda_(l) {}

  eos_t model_;
  std::vector<double> lambda_;
};

inline eos::eos_wrapper
init_eos(spec::config_py & config) {
  /*--------------------------------------------------------------------------*
     Equation of State
    *--------------------------------------------------------------------------*/
  const std::string name(config["eos"].cast<std::string>());

  if(name == "ideal")
    return singularity::IdealGas(
      config["gamma"].cast<double>() - 1, 2.0 /* FIXME */);
  else if(name == "spiner")
    return singularity::SpinerEOSDependsRhoSie(
      config["spiner_file"].cast<std::string>(),
      config["spiner_matid"].cast<std::string>());
  else if(name == "gruneisen")
    return singularity::Gruneisen(config["gruneisen_c0"].cast<double>(),
      config["gruneisen_s1"].cast<double>(),
      0., // s2
      0., // s3
      config["gruneisen_G0"].cast<double>(),
      0., // b
      0., // rho0
      config["gruneisen_T0"].cast<double>(),
      0., // P0
      config["gruneisen_Cv"].cast<double>(),
      0. // rho_max
    );
  else if(name == "jwl")
    return singularity::JWL(config["JWL_A"].cast<double>(),
      config["JWL_B"].cast<double>(),
      config["JWL_R1"].cast<double>(),
      config["JWL_R2"].cast<double>(),
      config["JWL_w"].cast<double>(),
      config["JWL_rho0"].cast<double>(),
      config["JWL_Cv"].cast<double>());

  flog_fatal("unsupported EOS" << name);
  return eos::eos_wrapper{};
}
inline void
init_eos_mm(std::vector<eos::eos_wrapper> & eos_mm, spec::config_py & config) {
  if(!config.contains("eos_list")) {
    flog_fatal("missing 'eos_list' section for multi-material EOS");
  }

  auto mats = config["eos_list"];

  if(!mats.is_sequence() || mats.empty()) {
    flog_fatal("'eos_list' must be a non-empty sequence");
  }

  eos_mm.clear();
  eos_mm.reserve(mats.size());

  for(std::size_t m = 0; m < mats.size(); ++m) {
    auto mat = mats[m];
    eos_mm.emplace_back(init_eos(mat));
  }
}

template<class Exec>
struct eos_dispatch {
  using exec_t = std::remove_cvref_t<Exec>;

  static constexpr bool on_device =
    std::is_same_v<exec_t, flecsi::exec::gpu> ||
    std::is_same_v<exec_t, flecsi::exec::accelerator>;

  using eos_type = eos_wrapper;

  std::array<eos_type, kMaxMaterials> models_{};
  std::size_t size_ = 0;
  bool initialized_ = false;

  explicit eos_dispatch(std::vector<eos_wrapper> const & eos_vec) {
    if(eos_vec.size() > kMaxMaterials) {
      throw std::runtime_error("too many materials for eos_dispatch");
    }

    size_ = eos_vec.size();

    for(std::size_t i = 0; i < size_; ++i) {
      if constexpr(on_device) {
        models_[i] = eos_wrapper::from_host_device(eos_vec[i]);
      }
      else {
        models_[i] = eos_vec[i];
      }
    }

    initialized_ = true;
  }

  ~eos_dispatch() {
    if constexpr(on_device) {
      if(initialized_) {
        for(std::size_t i = 0; i < size_; ++i) {
          models_[i].finalize_device();
        }
      }
    }

    initialized_ = false;
    size_ = 0;
  }

  FLECSI_INLINE_TARGET eos_type const & operator[](std::size_t mat) const {
    return models_[mat];
  }

  FLECSI_INLINE_TARGET std::size_t size() const {
    return size_;
  }
};

} // namespace eos

#endif // HYDRO_EOS_HH