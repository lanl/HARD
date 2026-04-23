#include "options.hh"
#include "state.hh"

#include "common/utils.hh"
#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/init.hh"
#include "hydro/tasks/initial_data/all_initial_data.hh"
#include "hydro/tasks/maxcharspeed.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/eos.hh"
#include "spec/tasks/boundaries/boundary.hh"
#include "spec/tasks/io.hh"
#include "spec/types.hh"

#include "../../actions/initialize/init.hh"

#include <spec/runtime.hh>

#include <flecsi/flog.hh>

namespace hard {

template<std::size_t D>
struct initialize {
  static void action(control_policy<state, D> & cp) {
    using namespace flecsi;
    using namespace common::utils;
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    spec::config_py config(opt::config.value());

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
      Mesh topology allocation
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

    if(config["problem"].cast<std::string>() == "sod") {
      execute<
        tasks::initial_data::shock<tasks::initial_data::shock_tubes::sod, D>>(
        flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.eos);
    }
    else if(config["problem"].cast<std::string>() == "rankine-hugoniot") {
      execute<tasks::initial_data::
          shock<tasks::initial_data::shock_tubes::rankine_hugoniot, D>>(
        flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.eos);
    }
    else if(config["problem"].cast<std::string>() == "leblanc") {
      execute<tasks::initial_data::
          shock<tasks::initial_data::shock_tubes::leblanc, D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.eos);
    }
    else if(config["problem"].cast<std::string>() == "acoustic-wave") {
      execute<tasks::initial_data::acoustic_wave<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.eos);
    }
    else if(config["problem"].cast<std::string>() == "kh-test") {
      execute<tasks::initial_data::kh_instability<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.eos);
    }
    // Rayleigh-Taylor setup
    else if(config["problem"].cast<std::string>() == "rt-test") {
      execute<tasks::initial_data::rt_instability<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.src_t.hydro.gravity_force(*s.m),
        s.icst.gravity_acc(*s.gt),
        s.cons.hydro.total_energy_density(*s.m),
        s.eos);
    }
    else if(config["problem"].cast<std::string>() == "sedov") {
      execute<tasks::initial_data::sedov_blast<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m));
    }
    else if(config["problem"].cast<std::string>() == "lw-implosion") {
      execute<tasks::initial_data::lw_implosion<D>>(flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        config["gamma"].cast<double>());
    }
    else {
      flog_fatal(
        "unsupported problem(" << config["problem"].cast<std::string>() << ")");
    } // if

    /*--------------------------------------------------------------------------*
      Initialize time advance.
      *--------------------------------------------------------------------------*/

    actions::init_timestep(s, sc, config);
  }
};

static const auto initialize_action =
  spec::register_action<control, state, initialize, cp::initialize>();

} // namespace hard
