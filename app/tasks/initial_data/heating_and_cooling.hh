
#pragma once

#include "../../constants.hh"
#include "../../options.hh"
#include "../../types.hh"
#include <cstddef>
#include <spec/utils.hh>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

//
// A test problem demonstrating fluid and radiation reaching to a thermal
// equilibrium.
//
// The fluid is uniform and initially at rest everywhere. Radiation energy
// density is also uniform everywhere.
//
// Ideally, since E_rad = const. and nonzero everywhere, this problem should run
// without any issues even if every radiation-hydro tasks are switched on.
//
template<std::size_t Dim>
auto
heating_and_cooling(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> particle_mass_a,
  const double gamma) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto const particle_mass = *particle_mass_a;

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const double mass_density_v =
    config["problem_parameters"]["fluid_mass_density"].as<double>();
  const double fluid_temperature =
    config["problem_parameters"]["fluid_temperature"].as<double>();
  const double radiation_temperature =
    config["problem_parameters"]["radiation_temperature"].as<double>();

  // Note : assuming ideal gas EOS
  const double fluid_internal_energy_density =
    mass_density_v * constants::cgs::boltzmann_constant * fluid_temperature /
    ((gamma - 1.0) * particle_mass);
  const double radiation_energy_density_v =
    constants::cgs::radiation_constant *
    spec::utils::sqr(spec::utils::sqr(radiation_temperature));

  if constexpr(Dim == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      mass_density(i) = mass_density_v;
      momentum_density(i).x = 0.0;

      total_energy_density(i) = fluid_internal_energy_density;
      radiation_energy_density(i) = radiation_energy_density_v;
    }
  }
  else {
    flog_fatal(
      "Radiative heating-cooling problem is only implemented for D = 1")
  }
}

} // namespace hard::tasks::initial_data
