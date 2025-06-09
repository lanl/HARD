
#pragma once

#include "../../types.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

namespace shock_tubes {

struct rankine_hugoniot {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.25;
  static constexpr double uR = 0.0;
  static constexpr double pR = 0.1795;

  static constexpr double x0 = 0.5;
}; // struct rankine_hugoniot

struct sod {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 0.0;
  static constexpr double pR = 0.1;

  static constexpr double x0 = 0.5;
}; // struct sod

struct leblanc {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 0.1;

  static constexpr double rR = 1.0e-3;
  static constexpr double uR = 0.0;
  static constexpr double pR = 1.0e-10;

  static constexpr double x0 = 0.5;
}; // struct leblanc

} // namespace shock_tubes

/*----------------------------------------------------------------------------*
  Shock Tube.
 *----------------------------------------------------------------------------*/

template<typename T, std::size_t D>
auto
shock(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<D>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  const eos::eos_wrapper & eos) {
  using hard::tasks::util::get_mdiota_policy;

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const auto x = m.template head<ax::x>(i);

      if(x < T::x0) {
        mass_density(i) = T::rL;
        momentum_density(i).x = T::rL * T::uL;
        const double e = util::find_sie(eos, T::rL, T::pL);
        total_energy_density(i) = T::rL * e + 0.5 * T::rL * (T::uL * T::uL);
      }
      else {
        mass_density(i) = T::rR;
        momentum_density(i).x = T::rR * T::uR;
        const double e = util::find_sie(eos, T::rR, T::pR);
        total_energy_density(i) = T::rR * e + 0.5 * T::rR * (T::uR * T::uR);
      } // if

      radiation_energy_density(i) = 0.0;

    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_yx) {
      auto [j, i] = ji;
      const auto x = m.template head<ax::x>(i);

      if(x < T::x0) {
        mass_density(i, j) = T::rL;
        momentum_density(i, j) = vec<D>(T::rL * T::uL);
        const double e = util::find_sie(eos, T::rL, T::pL);
        total_energy_density(i, j) = T::rL * e + 0.5 * T::rL * (T::uL * T::uL);
      }
      else {
        mass_density(i, j) = T::rR;
        momentum_density(i, j) = vec<D>(T::rR * T::uR);
        const double e = util::find_sie(eos, T::rR, T::pR);
        total_energy_density(i, j) = T::rR * e + 0.5 * T::rR * (T::uR * T::uR);
      } // if

      radiation_energy_density(i, j) = 0.0;
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_zyx) {
      auto [k, j, i] = kji;
      const auto x = m.template head<ax::x>(i);

      if(x < T::x0) {
        mass_density(i, j, k) = T::rL;
        momentum_density(i, j, k) = vec<D>(T::rL * T::uL);
        const double e = util::find_sie(eos, T::rL, T::pL);
        total_energy_density(i, j, k) =
          T::rL * e + 0.5 * T::rL * (T::uL * T::uL);
      }
      else {
        mass_density(i, j, k) = T::rR;
        momentum_density(i, j, k) = vec<D>(T::rR * T::uR);
        const double e = util::find_sie(eos, T::rR, T::pR);
        total_energy_density(i, j, k) =
          T::rR * e + 0.5 * T::rR * (T::uR * T::uR);
      } // if

      radiation_energy_density(i, j, k) = 0.0;
    }; // forall
  } // if
} // shock

} // namespace hard::tasks::initial_data
