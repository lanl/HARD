#ifndef HYDRO_EOS_HH
#define HYDRO_EOS_HH

#include "types.hh"

#include <flecsi/util/serialize.hh>

#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos.hpp>

#include <cmath>
#include <utility>

namespace eos {

struct eos_wrapper {
  using eos_t = singularity::Variant<singularity::IdealGas,
    singularity::SpinerEOSDependsRhoSie,
    singularity::Gruneisen,
    singularity::JWL>;

  eos_wrapper() {}
  template<typename M>
  eos_wrapper(M const & m) : model_(m), lambda_(model_.nlambda()) {}

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
  eos_wrapper(eos_t const & eos, std::vector<double> const & l)
    : model_(eos), lambda_(l) {}
  eos_t model_;
  std::vector<double> lambda_;
}; // struct eos_base

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
init_eos_mm(std::vector<eos::eos_wrapper> eos_mm, spec::config_py & config) {
  if(!config["eos_list"]) {
    flog_fatal("missing 'eos_list' section for multi-material EOS");
  }

  auto mats = config["eos_list"];

  if(!mats.IsSequence() || mats.size() == 0) {
    flog_fatal("'eos_list' must be a non-empty sequence");
  }

  eos_mm.clear();
  eos_mm.reserve(mats.size());

  for(std::size_t m = 0; m < mats.size(); ++m) {
    auto mat = mats[m];
    eos_mm.emplace_back(init_eos(mat));
  }
}

} // namespace eos

#endif // HYDRO_EOS_HH
