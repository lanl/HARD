#include "options.hh"
#include "rad/tasks/utils.hh"
#include "state.hh"
#include "types.hh"

#include "common/utils.hh"
#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/init.hh"
#include "hydro/tasks/maxcharspeed.hh"
#include "hydro/tasks/rhs.hh"
#include "rad/tasks/init.hh"
#include "rad/tasks/initial_data/all_initial_data.hh"
#include "spec/eos.hh"
#include "spec/tasks/boundaries/boundary.hh"
#include "spec/tasks/io.hh"

#include "../../actions/initialize/init.hh"

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

namespace hard {

template<std::size_t D>
void
initialize(control_policy<state, D> & cp) {
  using namespace flecsi;
  using namespace common::utils;
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  YAML::Node config = YAML::LoadFile(opt::config.value());

  /*--------------------------------------------------------------------------*
    Solver.
   *--------------------------------------------------------------------------*/

  s.rad.mgr.solver_settings.maxiter =
    config["linear_solver"]["maxiter"].IsDefined()
      ? config["linear_solver"]["maxiter"].as<int>()
      : 50;
  s.rad.mgr.solver_settings.rtol =
    config["linear_solver"]["rtol"].IsDefined()
      ? config["linear_solver"]["rtol"].as<double>()
      : 1e-12;
  s.rad.mgr.solver_settings.use_zero_guess =
    config["linear_solver"]["use_zero_guess"].IsDefined()
      ? config["linear_solver"]["use_zero_guess"].as<bool>()
      : true;
  s.rad.mgr.flecsolve_coarse_grid =
    config["linear_solver"]["flecsolve_coarse_grid"].IsDefined()
      ? config["linear_solver"]["flecsolve_coarse_grid"].as<bool>()
      : true;
  s.rad.mgr.jacobi_iterations =
    config["linear_solver"]["jacobi_iterations"].IsDefined()
      ? config["linear_solver"]["jacobi_iterations"].as<double>()
      : 100;

  s.rad.mgr.full_multigrid =
    config["linear_solver"]["full_multigrid"].IsDefined()
      ? config["linear_solver"]["full_multigrid"].as<bool>()
      : false;

  /*--------------------------------------------------------------------------*
    Global and color topology allocations.
   *--------------------------------------------------------------------------*/

  std::vector<double> time;
  std::vector<double> temperature;
  actions::init_topologies(s, sc, time, temperature, config);

  /*--------------------------------------------------------------------------*
    Set boundaries.
   *--------------------------------------------------------------------------*/

  auto bf = actions::init_boundaries(s, sc, time, temperature, config);

  /*--------------------------------------------------------------------------*
    Kappa.
    *--------------------------------------------------------------------------*/

  execute<tasks::init::kappa>(
    s.rad.icst.kappa(*s.gt), config["kappa"].as<double>());

  /*--------------------------------------------------------------------------*
    Adaptive FLD Check, Closure ID and Limiter ID
   *--------------------------------------------------------------------------*/

  // Default is limiter = 1 and closure = 3
  std::size_t ci = config["closure_id"].IsDefined()
                     ? config["closure_id"].as<std::size_t>()
                     : 3;
  std::size_t li = config["limiter_id"].IsDefined()
                     ? config["limiter_id"].as<std::size_t>()
                     : 1;
  sc.execute<tasks::init::closure_id>(s.rad.icst.closure_id(*s.gt), ci);
  sc.execute<tasks::init::limiter_id>(s.rad.icst.limiter_id(*s.gt), li);

  /*--------------------------------------------------------------------------*
    Gravity Acceleration
   *--------------------------------------------------------------------------*/

  vec<D> g(0.0);
  if(config["gravity_acc"].IsDefined()) {
    g[0] = config["gravity_acc"][0].as<double>();
    if constexpr(D > 1)
      g[1] = config["gravity_acc"][1].as<double>();
    if constexpr(D > 2)
      g[2] = config["gravity_acc"][2].as<double>();
  }
  execute<tasks::init::initialize_gravity_acc<D>>(s.icst.gravity_acc(*s.gt), g);

  /*--------------------------------------------------------------------------*
    Particle mass
   *--------------------------------------------------------------------------*/
  execute<tasks::init::particle_mass>(
    s.icst.particle_mass(*s.gt), config["mean_molecular_weight"].as<double>());

  /*--------------------------------------------------------------------------*
    Mesh topology allocation.
   *--------------------------------------------------------------------------*/

  actions::init_mesh(s, sc, bf, config);

  /*--------------------------------------------------------------------------*
    Equation of State
   *--------------------------------------------------------------------------*/

  actions::init_eos(s, sc, config);

  /*--------------------------------------------------------------------------*
    Initialize problem state.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::initialize_gravity_force<D>>(
    flecsi::exec::on, s.src_t.hydro.gravity_force(*s.m));

  // Ritchmyer-Meshkov works with both radiation on and off
  if(config["problem"].as<std::string>() == "richtmyer-meshkov") {
    execute<tasks::initial_data::richtmyer_meshkov<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "heating_and_cooling") {
    if(config["eos"].as<std::string>() != "ideal")
      flog_fatal("Heating and cooling test only supports Ideal Gas eos");
    execute<tasks::initial_data::heating_and_cooling<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.icst.particle_mass(*s.gt),
      config["gamma"].as<double>());
  }
  // Heating and Cooling for AFLD
  else if(config["problem"].as<std::string>() == "heating-cooling-afld") {
    if(config["eos"].as<std::string>() != "ideal")
      flog_fatal("Heating and cooling test only supports Ideal Gas eos");
    sc.execute<tasks::initial_data::heating_and_cooling_afld<D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.icst.particle_mass(*s.gt),
      config["gamma"].as<double>());
  }
  else if(config["problem"].as<std::string>() == "implosion") {
    execute<tasks::initial_data::implosion_forced_T<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.icst.temperature_boundary(*s.dense_topology),
      s.icst.particle_mass(*s.gt),
      config["gamma"].as<double>());
  }
  // FIXME: This problem has not been tested for correctness
  else if(config["problem"].as<std::string>() == "rad-rh") {
    execute<tasks::initial_data::
        rad_RH<tasks::initial_data::rad_shock::rad_rankine_hugoniot, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      config["gamma"].as<double>(),
      s.icst.particle_mass(*s.gt));
  }
  // Kelvin Helmholtz with radiation setup
  else if(config["problem"].as<std::string>() == "kh-rad-test") {

    execute<tasks::initial_data::kh_instability_rad<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.eos);
  }
  else {
    flog_fatal(
      "unsupported problem(" << config["problem"].as<std::string>() << ")");
  } // if

  /*--------------------------------------------------------------------------*
    Initialize time advance.
   *--------------------------------------------------------------------------*/

  actions::init_timestep(
    s, sc, config, std::vector{s.rad.cons.radiation_energy_density(*s.m)});

} // initialize

inline control<state, 1>::action<initialize<1>, cp::initialize> init_1d;
inline control<state, 2>::action<initialize<2>, cp::initialize> init_2d;
inline control<state, 3>::action<initialize<3>, cp::initialize> init_3d;

} // namespace hard
