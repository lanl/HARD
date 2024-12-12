
#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>
#include <yaml-cpp/yaml.h>

namespace flastro::tasks::initial_data {

//
template<std::size_t Dim>
auto
shadow_2d(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // adiabatic index, assuming ideal gas EOS
  auto const gamma = *gamma_a;
  const double mult = 1.0 / (gamma - 1.0);

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());
  const double rho_high =
    config["problem_parameters"]["mass_density_high"].as<double>();
  const double rho_low =
    config["problem_parameters"]["mass_density_low"].as<double>();
  const double pressure = config["problem_parameters"]["pressure"].as<double>();
  const double Erad_left =
    config["problem_parameters"]["radiation_energy_density_left"].as<double>();
  const double Erad_right =
    config["problem_parameters"]["radiation_energy_density_right"].as<double>();

  const double x0 = config["problem_parameters"]["x0"].as<double>();
  const double y0 = config["problem_parameters"]["y0"].as<double>();

  const double r0_squared = spec::utils::sqr(
    config["problem_parameters"]["cylinder_radius"].as<double>());

  if constexpr(Dim == 3) {

    forall(k, (m.template cells<ax::z, dm::quantities>()), "init_shadow_2d") {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {

          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);

          const double r_squared = (x - x0) * (x - x0) + (y - y0) * (y - y0);

          mass_density(i, j, k) = (r_squared < r0_squared) ? rho_high : rho_low;
          momentum_density(i, j, k).x = 0.0;
          momentum_density(i, j, k).y = 0.0;
          momentum_density(i, j, k).z = 0.0;
          total_energy_density(i, j, k) = pressure * mult;

          // radiation_energy_density(i, j, k) =
          // (1.0 - x) * Erad_left + x * Erad_right;
          radiation_energy_density(i, j, k) = 1e-5;

        } // for
      }
    };
  }
  else {
    flog_fatal("")
  }
}

} // namespace flastro::tasks::initial_data
