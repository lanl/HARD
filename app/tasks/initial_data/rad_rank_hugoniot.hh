#pragma once

#include "../../constants.hh"
#include "../../numerical_algorithms/root_finder.hh"
#include "../../types.hh"
#include <cmath>
#include <cstddef>
#include <spec/utils.hh>

namespace flastro::tasks::initial_data {

/*----------------------------------------------------------------------------*
  Radiation modified Rankine Hugoniot.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
rad_RH(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  typename field<vec<D>>::template accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  field<double>::accessor<rw, ro> Erad_a,
  single<double>::accessor<ro> gamma_a,
  single<double>::accessor<ro> particle_mass_a) {
  auto r = m.template mdcolex<is::cells>(r_a);
  auto ru = m.template mdcolex<is::cells>(ru_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);
  auto Erad = m.template mdcolex<is::cells>(Erad_a);
  auto const gamma = *gamma_a;
  auto const particle_mass = *particle_mass_a;
  const double mult = 1.0 / (gamma - 1.0);

  const double kb = flastro::constants::cgs::boltzmann_constant;
  const double a = flastro::constants::cgs::radiation_constant;

  using spec::utils::sqr;

  // Left state
  const double mass_density_left = 1.0e-2;
  const double velocity_left = 1.0e9;
  const double temperature_left = 1.0e4;
  const double pressure_left =
    (mass_density_left / particle_mass) * kb * temperature_left;
  const double internal_energy_left = mult * pressure_left;
  const double total_energy_left =
    internal_energy_left + 0.5 * mass_density_left * sqr(velocity_left);
  const double radiation_energy_left = a * sqr(sqr(temperature_left));

  // Right state: solve with an exact jump condition
  const double mass_density_right = 6.85847e-2;
  const double velocity_right =
    mass_density_left * velocity_left / mass_density_right;

  const double c4 = a / 3.0;
  const double c1 = mass_density_right * kb / particle_mass;
  const double c0 = mass_density_right * sqr(velocity_right) -
                    (pressure_left + mass_density_left * sqr(velocity_left) +
                      radiation_energy_left / 3.0);
  const double temperature_right =
    numerical_algorithms::root_finder::halleys_get_root(
      c4, c1, c0, 1e4 * temperature_left);
  const double pressure_right =
    (mass_density_right / particle_mass) * kb * temperature_right;
  const double internal_energy_right = mult * pressure_right;
  const double total_energy_right =
    internal_energy_right + 0.5 * mass_density_right * sqr(velocity_right);
  const double radiation_energy_right = a * sqr(sqr(temperature_right));

  flog(info) << " -- Radiative shock problem -- " << std::endl;
  flog(info) << "  e_left = " << total_energy_left
             << "  E_left = " << radiation_energy_left << "\n"
             << "  v_right = " << velocity_right
             << "  e_right = " << total_energy_right
             << "  T_right = " << temperature_right
             << "  E_right = " << radiation_energy_right << std::endl;

  // Location of the shockfront
  const double x0 = 5.0e4;

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const auto x = m.template center<ax::x>(i);

      if(x < x0) {
        r(i) = mass_density_left;
        ru(i).x = mass_density_left * velocity_left;
        rE(i) = total_energy_left;
        Erad(i) = radiation_energy_left;
      }
      else {
        r(i) = mass_density_right;
        ru(i).x = mass_density_right * velocity_right;
        rE(i) = total_energy_right;
        Erad(i) = radiation_energy_right;
      } // if
      if(rE(i) <= 0) {
        std::cout << "TotalE is negative for i = " << i << std::endl;
      }
    } // for
  }
  // else if constexpr(D == 2) {
  //   for(auto j : m.template cells<ax::y, dm::quantities>()) {
  //     for(auto i : m.template cells<ax::x, dm::quantities>()) {
  //       const auto x = m.template center<ax::x>(i);

  //       if(x < T::x0) {
  //         r(i, j) = T::rL;
  //         ru(i, j).x = T::rL * T::uL;
  //         ru(i, j).y = T::rL * T::vL;
  //         // rE(i, j) = T::EL;
  //         // Erad(i, j) = T::EradL;
  //         rE(i, j) = mult * kb * T::TL * T::rL / (particle_mass) +
  //                    (0.5 * T::rL * T::uL * T::uL);
  //         Erad(i, j) = a * T::TL * T::TL * T::TL * T::TL;
  //       }
  //       else {
  //         r(i, j) = T::rR;
  //         ru(i, j).x = T::rR * T::uR;
  //         ru(i, j).y = T::rR * T::vR;
  //         // rE(i, j) = T::ER;
  //         // Erad(i, j) = T::EradR;
  //         rE(i, j) = mult * kb * T::TR * T::rR / (particle_mass) +
  //                    (0.5 * T::rR * T::uR * T::uR);
  //         Erad(i, j) = a * T::TR * T::TR * T::TR * T::TR;
  //       } // if
  //     } // for
  //   } // for
  // }
  // else /* D == 3 */ {
  //   for(auto k : m.template cells<ax::z, dm::quantities>()) {
  //     for(auto j : m.template cells<ax::y, dm::quantities>()) {
  //       for(auto i : m.template cells<ax::x, dm::quantities>()) {
  //         const auto x = m.template center<ax::x>(i);

  //         if(x < T::x0) {
  //           r(i, j, k) = T::rL;
  //           ru(i, j, k).x = T::rL * T::uL;
  //           ru(i, j, k).y = T::rL * T::vL;
  //           ru(i, j, k).z = T::rL * T::wL;
  //           // rE(i, j, k) = T::EL;
  //           // Erad(i, j, k) = T::EradL;
  //           rE(i, j, k) = mult * kb * T::TL * T::rL / (particle_mass) +
  //                         (0.5 * T::rL * T::uL * T::uL);
  //           Erad(i, j, k) = a * T::TL * T::TL * T::TL * T::TL;
  //         }
  //         else {
  //           r(i, j, k) = T::rR;
  //           ru(i, j, k).x = T::rR * T::uR;
  //           ru(i, j, k).y = T::rR * T::vR;
  //           ru(i, j, k).z = T::rR * T::wR;
  //           // rE(i, j, k) = T::ER;
  //           // Erad(i, j, k) = T::EradR;
  //           rE(i, j, k) = mult * kb * T::TR * T::rR / (particle_mass) +
  //                         (0.5 * T::rR * T::uR * T::uR);
  //           Erad(i, j, k) = a * T::TR * T::TR * T::TR * T::TR;
  //         } // if
  //       } // for
  //     } // for
  //   } // for
  // } // if

} // rad_RH

} // namespace flastro::tasks::initial_data
