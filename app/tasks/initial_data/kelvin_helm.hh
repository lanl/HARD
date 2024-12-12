
#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>

namespace flastro::tasks::initial_data {

/*----------------------------------------------------------------------------*
  Kelvin-Helmholtz Instability.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
kh_instability(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a) {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto const gamma = *gamma_a;
  const double mult = 1.0 / (gamma - 1.0);

  using spec::utils::sqr;

  // density, velocity and pressure for two strips of fluids
  const double rL = 1.0;
  const double uL = -0.5;
  const double pL = 2.5;

  const double rH = 2.0;
  const double uH = 0.5;
  const double pH = 2.5;

  // Perturbation
  const double velocity_perturbation_amplitude = 0.1;
  const double N = 2;
  const double wavenumber = N * 2.0 * M_PI;
  const double perturb_width = 0.05 / sqrt(2.0);
  const double gaussian_factor = 0.5 / sqr(perturb_width);

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    forall(j, (m.template cells<ax::y, dm::quantities>()), "init_kh_2d") {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        const double v_y = velocity_perturbation_amplitude *
                           sin(wavenumber * x) *
                           (exp(-gaussian_factor * sqr(y - 0.25)) +
                             exp(-gaussian_factor * sqr(y - 0.75)));

        // initialize two different density and velocity fluids
        if(std::abs(y - 0.5) > 0.25) {
          mass_density(i, j) = rL;
          momentum_density(i, j).x = rL * uL;
          momentum_density(i, j).y = rL * v_y;
          total_energy_density(i, j) =
            mult * pL + 0.5 * rL * (sqr(uL) + sqr(v_y));
        }
        else {
          mass_density(i, j) = rH;
          momentum_density(i, j).x = rH * uH;
          momentum_density(i, j).y = rH * v_y;
          total_energy_density(i, j) =
            mult * pH + 0.5 * rH * (sqr(uH) + sqr(v_y));
        } // if

        radiation_energy_density(i, j) = 0.0;

      } // for
    }; // forall
  }
  else /* D == 3 */ {
    forall(k, (m.template cells<ax::z, dm::quantities>()), "init_kh_3d") {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {

          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);

          const double v_y = velocity_perturbation_amplitude *
                             sin(wavenumber * x) *
                             (exp(-gaussian_factor * sqr(y - 0.25)) +
                               exp(-gaussian_factor * sqr(y - 0.75)));

          // initialize two different density and velocity fluids
          if(std::abs(y - 0.5) > 0.25) {
            mass_density(i, j, k) = rL;
            momentum_density(i, j, k).x = rL * uL;
            momentum_density(i, j, k).y = rL * v_y;
            total_energy_density(i, j, k) =
              mult * pL + 0.5 * rL * (sqr(uL) + sqr(v_y));
          }
          else {
            mass_density(i, j, k) = rH;
            momentum_density(i, j, k).x = rH * uH;
            momentum_density(i, j, k).y = rH * v_y;
            total_energy_density(i, j, k) =
              mult * pH + 0.5 * rH * (sqr(uH) + sqr(v_y));
          } // if

          momentum_density(i, j, k).z = 0.0;
          radiation_energy_density(i, j, k) = 0.0;
        }
      }
    };
  }
} //  kh_instability

} // namespace flastro::tasks::initial_data
