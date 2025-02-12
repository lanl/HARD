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
    singularity::Gruneisen>;

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

} // namespace eos

template<>
struct flecsi::util::serial::traits<eos::eos_wrapper> {
  using type = eos::eos_wrapper;
  template<class P>
  static void put(P & p, type const & cw) {
    auto & w = const_cast<type &>(cw);
    auto [n, b] = w.model_.Serialize();
    serial::put(p, n);
    mempcpy(p, b, n);
    free(b);
    serial::put(p, w.lambda_);
  }
  static type get(std::byte const *& p) {
    const auto n = serial::get<std::size_t>(p);
    typename eos::eos_wrapper::eos_t model;
    model.DeSerialize(const_cast<char *>(reinterpret_cast<char const *>(p)));
    p += n;
    return {model, serial::get<std::vector<double>>(p)};
  }
};

#endif // HYDRO_EOS_HH
