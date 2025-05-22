
#pragma once

#include "../../options.hh"
#include "../../types.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>
#include <cstring>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

//
// An acoustic wave set-up
//
template<std::size_t Dim>
auto
acoustic_wave(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  const eos::eos_wrapper & eos) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  YAML::Node config = YAML::LoadFile(opt::config.value());

  // Problem parameters
  // Equilibrium values
  const double r0{
    config["problem_parameters"]["r0"].as<double>()}; // Equilibrium density
  const double p0{
    config["problem_parameters"]["p0"].as<double>()}; // Equilibrium pressure

  // Perturbation amplitudes
  const double rA{
    config["problem_parameters"]["amplitude"].as<double>()}; // Density
  const double uA{
    config["problem_parameters"]["amplitude"].as<double>()}; // Velocity

  // Sound speed
  const double cs{sqrt(config["gamma"].as<double>() * p0 / r0)};

  // Define the inverse of the wave number
  const double scale{config["problem_parameters"]["scale"].as<double>()};

  // sine_quad is the volume average of the sine per cell,
  // for sine wave fvm initialization
  const double k{2 * M_PI * scale};
  auto sine_quad = [k](double x0, double x1) {
    return (cos(k * x0) - cos(k * x1)) / (k * (x1 - x0));
  };

  //
  // Only 1D version has been implemented.
  //
  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "init_acoustic_1d") {
      const auto x0{m.template head<ax::x>(i)};
      const auto x1{m.template tail<ax::x>(i)};
      const auto x{m.template center<ax::x>(i)};
      const double ux{cs * uA * sine_quad(x0, x1)};
      mass_density(i) = r0 * (1 + rA * sine_quad(x0, x1));

      momentum_density(i).x = mass_density(i) * ux;
      const double e = util::find_sie(eos, mass_density(i), p0);
      total_energy_density(i) =
        mass_density(i) * e + 0.5 * mass_density(i) * utils::sqr(ux);

      radiation_energy_density(i) = 0.0;

    }; // forall
  }
  else {
    flog_fatal("Acoustic wave problem for D >= 2 is not implemented")
  }
} // acoustic_wave

} // namespace hard::tasks::initial_data
