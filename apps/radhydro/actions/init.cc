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
#include "spec/types.hh"

#include "../../actions/initialize/init.hh"

#include <flecsi/flog.hh>
#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct initialize {

  static void action(control_policy<state, D> & cp) {
    using namespace flecsi;
    using namespace common::utils;
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    /*--------------------------------------------------------------------------*
      Solver.
     *--------------------------------------------------------------------------*/

    spec::config_py config(opt::config.value());

    s.rad.mgr.solver_settings.maxiter =
      config["linear_solver"].contains("maxiter")
        ? config["linear_solver"]["maxiter"].cast<int>()
        : 50;
    s.rad.mgr.solver_settings.rtol =
      config["linear_solver"].contains("rtol")
        ? config["linear_solver"]["rtol"].cast<double>()
        : 1e-12;
    s.rad.mgr.solver_settings.use_zero_guess =
      config["linear_solver"].contains("use_zero_guess")
        ? config["linear_solver"]["use_zero_guess"].cast<bool>()
        : true;
    s.rad.mgr.flecsolve_coarse_grid =
      config["linear_solver"].contains("flecsolve_coarse_grid")
        ? config["linear_solver"]["flecsolve_coarse_grid"].cast<bool>()
        : true;
    s.rad.mgr.jacobi_iterations =
      config["linear_solver"].contains("jacobi_iterations")
        ? config["linear_solver"]["jacobi_iterations"].cast<double>()
        : 100;
    s.rad.mgr.full_multigrid =
      config["linear_solver"].contains("full_multigrid")
        ? config["linear_solver"]["full_multigrid"].cast<bool>()
        : false;

    /*--------------------------------------------------------------------------*
      Global and color topology allocations.
     *--------------------------------------------------------------------------*/

    std::vector<double> time;
    std::vector<double> temperature;
    actions::init_topologies(s, sc, time, temperature);

    /*--------------------------------------------------------------------------*
      Set boundaries.
     *--------------------------------------------------------------------------*/

    auto bf = actions::init_boundaries(s, sc, time, temperature, config);

    /*--------------------------------------------------------------------------*
      Kappa.
      *--------------------------------------------------------------------------*/
    execute<tasks::init::kappa>(
      s.rad.icst.kappa(*s.gt), config["kappa"].cast<double>());

    /*--------------------------------------------------------------------------*
      Adaptive FLD Check, Closure ID and Limiter ID
     *--------------------------------------------------------------------------*/

    // Default is limiter = 1 and closure = 3
    const std::size_t ci = config.contains("closure_id")
                             ? config["closure_id"].cast<std::size_t>()
                             : 3;
    const std::size_t li = config.contains("limiter_id")
                             ? config["limiter_id"].cast<std::size_t>()
                             : 1;
    sc.execute<tasks::init::closure_id>(s.rad.icst.closure_id(*s.gt), ci);
    sc.execute<tasks::init::limiter_id>(s.rad.icst.limiter_id(*s.gt), li);

    /*--------------------------------------------------------------------------*
      Gravity Acceleration
     *--------------------------------------------------------------------------*/

    vec<D> g(0.0);
    if(config.contains("gravity_acc")) {
      g[0] = config["gravity_acc"][0].cast<double>();
      if constexpr(D > 1)
        g[1] = config["gravity_acc"][1].cast<double>();
      if constexpr(D > 2)
        g[2] = config["gravity_acc"][2].cast<double>();
    }
    execute<tasks::init::initialize_gravity_acc<D>>(
      s.icst.gravity_acc(*s.gt), g);

    /*--------------------------------------------------------------------------*
      Particle mass
     *--------------------------------------------------------------------------*/
    execute<tasks::init::particle_mass>(s.icst.particle_mass(*s.gt),
      config["mean_molecular_weight"].cast<double>());

    /*--------------------------------------------------------------------------*
      Mesh topology allocation.
     *--------------------------------------------------------------------------*/

    actions::init_mesh(s, sc, bf, config);

    /*--------------------------------------------------------------------------*
      Equation of State
     *--------------------------------------------------------------------------*/

    actions::init_eos(s, config);

    /*--------------------------------------------------------------------------*
      Initialize problem state.
     *--------------------------------------------------------------------------*/

    execute<tasks::init::initialize_gravity_force<D>>(
      flecsi::exec::on, s.src_t.hydro.gravity_force(*s.m));

    // Ritchmyer-Meshkov works with both radiation on and off
    if(config["problem"].cast<std::string>() == "richtmyer-meshkov") {
      execute<tasks::initial_data::richtmyer_meshkov<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.rad.cons.radiation_energy_density(*s.m),
        s.eos);
    }
    else if(config["problem"].cast<std::string>() == "heating_and_cooling") {
      if(config["eos"].cast<std::string>() != "ideal")
        flog_fatal("Heating and cooling test only supports Ideal Gas eos");
      execute<tasks::initial_data::heating_and_cooling<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.rad.cons.radiation_energy_density(*s.m),
        s.icst.particle_mass(*s.gt),
        config["gamma"].cast<double>());
    }
    // Heating and Cooling for AFLD
    else if(config["problem"].cast<std::string>() == "heating-cooling-afld") {
      if(config["eos"].cast<std::string>() != "ideal")
        flog_fatal("Heating and cooling test only supports Ideal Gas eos");
      sc.execute<tasks::initial_data::heating_and_cooling_afld<D>>(
        flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.rad.cons.radiation_energy_density(*s.m),
        s.icst.particle_mass(*s.gt),
        config["gamma"].cast<double>());
    }
    else if(config["problem"].cast<std::string>() == "implosion") {
      execute<tasks::initial_data::implosion_forced_T<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.rad.cons.radiation_energy_density(*s.m),
        s.icst.temperature_boundary(*s.dense_topology),
        s.icst.particle_mass(*s.gt),
        config["gamma"].cast<double>());
    }
    // FIXME: This problem has not been tested for correctness
    else if(config["problem"].cast<std::string>() == "rad-rh") {
      execute<tasks::initial_data::
          rad_RH<tasks::initial_data::rad_shock::rad_rankine_hugoniot, D>>(
        flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.rad.cons.radiation_energy_density(*s.m),
        config["gamma"].cast<double>(),
        s.icst.particle_mass(*s.gt));
    }
    // Kelvin Helmholtz with radiation setup
    else if(config["problem"].cast<std::string>() == "kh-rad-test") {

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
        "unsupported problem(" << config["problem"].cast<std::string>() << ")");
    } // if

    /*--------------------------------------------------------------------------*
      Initialize time advance.
     *--------------------------------------------------------------------------*/

    actions::init_timestep(
      s, sc, config, std::vector{s.rad.cons.radiation_energy_density(*s.m)});
  }
};

static const auto initialize_action =
  spec::register_action<control, state, initialize, cp::initialize>();

} // namespace hard
