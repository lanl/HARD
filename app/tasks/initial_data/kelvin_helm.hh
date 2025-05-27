
#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

/*----------------------------------------------------------------------------*
  Kelvin-Helmholtz Instability.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
kh_instability(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  const eos::eos_wrapper & eos) {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  const double wavenumber = 2.0 * M_PI;
  const double N = 4;

  // density, velocity and pressure for two fluids
  const double rL = 1.0;
  const double uL = -0.5;
  const double vL = 0.0;
  const double pL = 2.5;

  const double rH = 2.0;
  const double uH = 0.5;
  const double vH = 0.0;
  const double pH = 2.5;

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    s.executor().forall(j, (m.template cells<ax::x, dm::quantities>())) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        // initialize two different density and velocity fluids
        if(std::abs(y - 0.5) > 0.25) {
          mass_density(i, j) = rL;
          momentum_density(i, j).x = rL * uL;
          momentum_density(i, j).y = rL * vL;
          const double e = util::find_sie(eos, rL, pL);
          total_energy_density(i, j) = rL * e + 0.5 * rL * (vL * vL);
        }
        else {
          mass_density(i, j) = rH;
          momentum_density(i, j).x = rH * uH;
          momentum_density(i, j).y = rH * vH;
          const double e = util::find_sie(eos, rH, pH);
          total_energy_density(i, j) = rH * e + 0.5 * rH * (vH * vH);
        } // if

        radiation_energy_density(i, j) = 0.0;

        // velocity perturbations in the Y-direction
        if(std::abs(y - 0.25) < 0.1) {
          momentum_density(i, j).y = 0.05 * sin(N * wavenumber * x);
        }
        if(std::abs(y - 0.75) < 0.1) {
          momentum_density(i, j).y = 0.05 * sin(N * wavenumber * x);
        }
      } // for
    }; // forall
  }
  else /* D == 3 */ {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 3 is not implemented")
  } // if
} //  kh_instability

} // namespace hard::tasks::initial_data
