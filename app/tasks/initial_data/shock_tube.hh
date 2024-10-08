
#pragma once

#include "../../types.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

namespace shock_tubes {

struct rankine_hugoniot {
  static constexpr double rL = 2.299156e-01;
  static constexpr double uL = 1.270775e+00;
  static constexpr double vL = 0.0;
  static constexpr double wL = 0.0;
  static constexpr double pL = 6.700238e-01;

  static constexpr double rR = 1.250000e-01;
  static constexpr double uR = 0.000000e+00;
  static constexpr double vR = 0.0;
  static constexpr double wR = 0.0;
  static constexpr double pR = 2.276638e-01;

  static constexpr double x0 = 0.5;
  static constexpr double y0 = 0.5;
}; // struct rankine_hugoniot

struct sod {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double vL = 0.0;
  static constexpr double wL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 0.0;
  static constexpr double vR = 0.0;
  static constexpr double wR = 0.0;
  static constexpr double pR = 0.1;

  static constexpr double x0 = 0.5;
  static constexpr double y0 = 0.5;
}; // struct sod

} // namespace shock_tubes

/*----------------------------------------------------------------------------*
  Shock Tube.
 *----------------------------------------------------------------------------*/

template<typename T, std::size_t D>
auto
shock(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<D>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a) {
  using hard::tasks::util::get_mdiota_policy;

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "init_1d") {
      const double gamma = *gamma_a;
      const double mult = 1.0 / (gamma - 1.0);

      const auto x = m.template head<ax::x>(i);

      if(x < T::x0) {
        mass_density(i) = T::rL;
        momentum_density(i).x = T::rL * T::uL;
        total_energy_density(i) =
          mult * T::pL + 0.5 * T::rL * (utils::sqr(T::uL));
      }
      else {
        mass_density(i) = T::rR;
        momentum_density(i).x = T::rR * T::uR;
        total_energy_density(i) =
          mult * T::pR + 0.5 * T::rR * (utils::sqr(T::uR));
      } // if

      radiation_energy_density(i) = 0.0;

    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_yx, "init_shock") {
      auto [j, i] = ji;
      const double gamma = *gamma_a;
      const double mult = 1.0 / (gamma - 1.0);

      const auto x = m.template head<ax::x>(i);
      const auto y = m.template head<ax::y>(j);

      double radius =
        sqrt((x - T::x0) * (x - T::x0) + (y - T::y0) * (y - T::y0));

      // if(x < T::x0) {
      if(radius < 0.2) {
        mass_density(i, j) = T::rL;
        momentum_density(i, j).x = T::rL * T::uL;
        momentum_density(i, j).y = T::rL * T::vL;
        total_energy_density(i, j) =
          mult * T::pL + 0.5 * T::rL * (utils::sqr(T::uL) + utils::sqr(T::vL));
      }
      else {
        mass_density(i, j) = T::rR;
        momentum_density(i, j).x = T::rR * T::uR;
        momentum_density(i, j).y = T::rR * T::vR;
        total_energy_density(i, j) =
          mult * T::pR + 0.5 * T::rR * (utils::sqr(T::uR) + utils::sqr(T::vR));
      } // if

      radiation_energy_density(i, j) = 0.0;
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_zyx, "init_shock") {
      auto [k, j, i] = kji;
      const double gamma = *gamma_a;
      const double mult = 1.0 / (gamma - 1.0);

      const auto x = m.template head<ax::x>(i);

      if(x < T::x0) {
        mass_density(i, j, k) = T::rL;
        momentum_density(i, j, k).x = T::rL * T::uL;
        momentum_density(i, j, k).y = T::rL * T::vL;
        momentum_density(i, j, k).z = T::rL * T::wL;
        total_energy_density(i, j, k) =
          mult * T::pL +
          0.5 * T::rL *
            (utils::sqr(T::uL) + utils::sqr(T::vL) + utils::sqr(T::wL));
      }
      else {
        mass_density(i, j, k) = T::rR;
        momentum_density(i, j, k).x = T::rR * T::uR;
        momentum_density(i, j, k).y = T::rR * T::vR;
        momentum_density(i, j, k).z = T::rR * T::wR;
        total_energy_density(i, j, k) =
          mult * T::pR +
          0.5 * T::rR *
            (utils::sqr(T::uR) + utils::sqr(T::vR) + utils::sqr(T::wR));
      } // if

      radiation_energy_density(i, j, k) = 0.0;
    }; // forall
  } // if
} // shock

} // namespace hard::tasks::initial_data
