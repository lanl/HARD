
#pragma once

#include "../../constants.hh"
#include "../../types.hh"
#include <cstddef>
#include <spec/utils.hh>
#include <yaml-cpp/yaml.h>

namespace flastro::tasks::initial_data {

//
template<std::size_t Dim>
auto
rayleigh_taylor(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a,
  single<double>::accessor<ro> particle_mass_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto const gamma = *gamma_a;
  auto const particle_mass = *particle_mass_a;

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const double boundary_z =
    config["problem_parameters"]["boundary_z"].as<double>();

  const double density_upper =
    config["problem_parameters"]["density_upper"].as<double>();
  const double density_lower =
    config["problem_parameters"]["density_lower"].as<double>();
  const double p0 = config["problem_parameters"]["p0"].as<double>();
  const double velocity_amplitude =
    config["problem_parameters"]["velocity_amplitude"].as<double>();
  const double Erad =
    config["problem_parameters"]["radiation_energy_density"].as<double>();
  const double perturb_width =
    config["problem_parameters"]["pertub_width"].as<double>();
  const double wavelength =
    config["problem_parameters"]["wavelength"].as<double>();

  const double mult = 1.0 / (gamma - 1.0);

  const double g = 1.0;

  const double zmax = 20.0;

  const double gaussian_factor = 1.0 / (0.5 * perturb_width * perturb_width);

  constexpr size_t n_modes = 6;
  std::array<double, 2 * n_modes> random_modes;
  for(size_t i = 0; i < 2 * n_modes; i++) {
    random_modes.at(i) = -1.0 + 2.0 * static_cast<double>(rand()) /
                                  (static_cast<double>(RAND_MAX));
  }

  if constexpr(Dim == 3) {
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {

          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);
          const auto z = m.template center<ax::z>(k);

          if(z > boundary_z) {
            mass_density(i, j, k) = density_upper;
            const double pressure = p0 + density_upper * g * (zmax - z);
            total_energy_density(i, j, k) = mult * pressure;
          }
          else {
            mass_density(i, j, k) = density_lower;
            const double pressure =
              p0 + g * (density_upper * (zmax - boundary_z) +
                         density_lower * (boundary_z - z));
            total_energy_density(i, j, k) = mult * pressure;
          }

          double vz = 0.0;
          for(size_t n = 0; n < n_modes; ++n) {
            vz += velocity_amplitude * random_modes.at(2 * n) *
                  cos((n + 1) * 2.0 * M_PI * x / wavelength);
            vz += velocity_amplitude * random_modes.at(2 * n + 1) *
                  sin((n + 1) * 2.0 * M_PI * x / wavelength);
          }
          vz *= exp(-gaussian_factor * spec::utils::sqr(z - boundary_z));

          momentum_density(i, j, k).z = vz * mass_density(i, j, k);
          momentum_density(i, j, k).x = 0.0;
          momentum_density(i, j, k).y = 0.0;

          radiation_energy_density(i, j, k) = Erad;
        }
      }
    }
  }
  else {
    flog_fatal(
      "Rayleigh-Taylor instability problem is only implemented for D = 3")
  }
}

} // namespace flastro::tasks::initial_data
