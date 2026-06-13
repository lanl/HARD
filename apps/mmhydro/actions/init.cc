#include "options.hh"
#include "state.hh"

#include "common/utils.hh"
#include "gfm/tasks/cons2prim.hh"
#include "gfm/tasks/initial_data/all_initial_data.hh"
#include "gfm/tasks/ls_core.hh"
#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/init.hh"
#include "hydro/tasks/initial_data/all_initial_data.hh"
#include "hydro/tasks/maxcharspeed.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/eos.hh"
#include "spec/tasks/boundaries/boundary.hh"
#include "spec/tasks/io.hh"

#include "../../actions/initialize/init.hh"

#include <flecsi/flog.hh>

namespace hard {

template<std::size_t D>
void
initialize(control_policy<state, D> & cp) {
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
  execute<tasks::init::initialize_gravity_acc<D>>(s.icst.gravity_acc(*s.gt), g);

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

  eos::init_eos_mm(s.eos_mm, config);

  /*--------------------------------------------------------------------------*
    Initialize problem state.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::initialize_gravity_force<D>>(
    flecsi::exec::on, s.src_t.hydro.gravity_force(*s.m));

  auto problem = config["problem"].cast<std::string>();

  if(problem == "sod") {
    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::sod, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(problem == "sod_mm") {
    execute<tasks::initial_data::
        shock_mm<tasks::initial_data::shock_tubes::sod_mm, D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.ls.material(*s.m),
      s.ls.frozen_material(*s.m),
      s.ls.psi(*s.m),
      s.extp.q_scalar(*s.m),
      s.extp.q_scalar_known(*s.m),
      s.extp.q_ext_plus(*s.m),
      s.extp.q_ext_minus(*s.m),
      s.eos_mm);
  }
  else if(problem == "rankine-hugoniot") {
    execute<tasks::initial_data::
        shock<tasks::initial_data::shock_tubes::rankine_hugoniot, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(problem == "leblanc") {
    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::leblanc, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(problem == "acoustic-wave") {
    execute<tasks::initial_data::acoustic_wave<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(problem == "kh-test") {
    execute<tasks::initial_data::kh_instability<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  // Rayleigh-Taylor setup
  else if(problem == "rt-test") {
    execute<tasks::initial_data::rt_instability<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.src_t.hydro.gravity_force(*s.m),
      s.icst.gravity_acc(*s.gt),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(problem == "sedov") {
    execute<tasks::initial_data::sedov_blast<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m));
  }
  else if(problem == "lw-implosion") {
    execute<tasks::initial_data::lw_implosion<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      config["gamma"].cast<double>());
  }
  else {
    flog_fatal("unsupported problem(" << problem << ")");
  } // if

  /*--------------------------------------------------------------------------*
    Output method configuration.
    *--------------------------------------------------------------------------*/

  // Default output method is csv (set in state initialization)
  if(config.contains("output_method")) {
    auto method = config["output_method"].cast<std::string>();
    if(method == "vti") {
      s.output_method = spec::om::vti;
    }
#ifdef HARD_ENABLE_HDF5
    else if(method == "xdmf") {
      s.output_method = spec::om::xdmf;
    }
#endif
  }

  /*--------------------------------------------------------------------------*
    Initialize time advance.
    *--------------------------------------------------------------------------*/

  actions::init_timestep(s, sc, config);

} // initialize

inline control<state, 1>::action<initialize<1>, cp::initialize> init_1d;
inline control<state, 2>::action<initialize<2>, cp::initialize> init_2d;
inline control<state, 3>::action<initialize<3>, cp::initialize> init_3d;

} // namespace hard
