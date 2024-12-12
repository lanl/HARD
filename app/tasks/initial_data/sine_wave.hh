
#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>

namespace flastro::tasks::initial_data {

//
// A sinusoidal density profile being advected with a constant velocity.
// This test problem can be used to check the numerical convergence
//
template<std::size_t Dim>
auto
sine_wave(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // Problem parameters
  const double mass_density_amplitude = 0.2;
  const double wavenumber = 2.0 * M_PI; // wavelength = 1.0
  const double velocity_x =
    0.2; // later we might want to extend this to multi-D
  const double pressure = 1.0;

  //
  // Only 1D version has been implemented.
  //
  if constexpr(Dim == 1) {
    forall(
      i, (m.template cells<ax::x, dm::quantities>()), "init_sine_wave_1d") {
      auto const gamma = *gamma_a;
      const double mult = 1.0 / (gamma - 1.0);

      const auto x = m.template center<ax::x>(i);

      mass_density(i) = 1.0 + mass_density_amplitude * sin(wavenumber * x);
      momentum_density(i).x = mass_density(i) * velocity_x;
      total_energy_density(i) =
        mult * pressure + 0.5 * mass_density(i) * utils::sqr(velocity_x);

      radiation_energy_density(i) = 0.0;

    }; // forall
  }
  else {
    flog_fatal("Sine wave problem for D >= 2 is not implemented")
  }
}

} // namespace flastro::tasks::initial_data
