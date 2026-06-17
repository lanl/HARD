#ifndef HARD_MODULES_RAD_TASKS_INITIAL_IMPLOSION_FORCED_T_HH
#define HARD_MODULES_RAD_TASKS_INITIAL_IMPLOSION_FORCED_T_HH

#include "../../constants.hh"
#include "options.hh"
#include "types.hh"
#include <../modules/spec/utils.hh>

namespace hard::tasks::initial_data {

// A test problem that implodes inward with external temperature source
// from FDS input.
//
// This is only 1D problem for ICP-mockup type

struct implosion {
  static constexpr double ri = 1.0;
  static constexpr double pi = 1.0;
  static constexpr double ui = 1.0;
  static constexpr double Ti = 1.0;
}; // struct implosion

template<std::size_t Dim>
auto
implosion_forced_T(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  field<double>::accessor<ro> temperature_boundary_a,
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
  spec::config_py config(opt::config.value());
  const double mass_density_v =
    config["problem_parameters"]["fluid_mass_density"].cast<double>();
  const double fluid_temperature =
    config["problem_parameters"]["fluid_temperature"].cast<double>();

  // Constant radiation temperature in the domain
  const double radiation_temperature =
    config["problem_parameters"]["radiation_temperature"].cast<double>();

  // Note : assuming ideal gas EOS

  using constants::cgs::boltzmann_constant;
  using constants::cgs::radiation_constant;
  using spec::utils::sqr;

  const double fluid_internal_energy_density =
    mass_density_v * boltzmann_constant * fluid_temperature /
    ((gamma - 1.0) * particle_mass);

  const double radiation_energy_density_v =
    radiation_constant * sqr(sqr(radiation_temperature));

  if constexpr(Dim == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      mass_density(i) = mass_density_v;
      momentum_density(i).x() = 0.0;

      total_energy_density(i) = fluid_internal_energy_density;
      radiation_energy_density(i) = radiation_energy_density_v;
    };
  }
  else {
    flog_fatal("Implosion problem is only implemented for D == 1")
  }
}

} // namespace hard::tasks::initial_data

#endif // HARD_MODULES_RAD_TASKS_INITIAL_IMPLOSION_FORCED_T_HH
