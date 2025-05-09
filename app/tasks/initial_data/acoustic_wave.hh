
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

#ifdef ENABLE_RADIATION
  flog_fatal("Acoustic wave must be built with ENABLE_RADIATION=OFF");
#endif

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

  // Define the gaussian lambda
  const double sigma{config["problem_parameters"]["sigma"].as<double>()};
  const double x0{config["problem_parameters"]["x0"].as<double>()};
  auto gaussian = [sigma, x0](
                    double x) { return exp(-0.5 * pow((x - x0) / sigma, 2)); };

  const double k{2.0 * M_PI};
  auto sine = [k](double x) { return sin(k * x); };

  // NOTE: Allow for choice between gaussian and sine
  auto & initial_shape = gaussian;
  const std::string init{
    config["problem_parameters"]["init"].as<std::string>()};

  //
  // Only 1D version has been implemented.
  //
  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "init_acoustic_1d") {
      const auto x{m.template center<ax::x>(i)};
      const double ux{cs * uA * (init == "gaussian" ? gaussian(x) : sine(x))};

      mass_density(i) = r0 * (1 + rA * initial_shape(x));
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
