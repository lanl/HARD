
#pragma once

#include "../../types.hh"
#include <cstddef>
// #include <yaml-cpp/yaml.h>

namespace flastro::tasks::initial_data {

template<std::size_t Dim>
auto
su_olson(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());
  const double e_floor = config["problem_parameters"]["e_floor"].as<double>();
  const double E_floor = config["problem_parameters"]["E_floor"].as<double>();

  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "init_su_olson") {
      // auto const gamma = *gamma_a;
      // const double mult = 1.0 / (gamma - 1.0);

      mass_density(i) = 1.0;
      momentum_density(i).x = 0.0;
      total_energy_density(i) = e_floor;
      radiation_energy_density(i) = E_floor;
    };
  }
  else {
    flog_fatal("Wrong dimension for Su-Olson problem");
  }

  //
  // ....
  //
}

} // namespace flastro::tasks::initial_data
