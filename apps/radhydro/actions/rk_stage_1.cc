#include "rad/tasks/utils.hh"
#include "state.hh"

#include "rad/tasks/interface_fluxes.hh"
#include "rad/tasks/rad.hh"

#include "../../actions/rk_stages/rk_stage_1.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct hydro_rk_advance_1 {
  static void action(control_policy<state, D> & cp) {

    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    sc.execute<tasks::rad::explicit_source_update<D>>(flecsi::exec::on,
      *s.m,
      s.prim.velocity(*s.m),
      s.rad.src_t.radiation_force(*s.m),
      s.rad.src_t.radiation_pressure_tensor(*s.m),
      s.velocity_gradient(*s.m),
      //
      s.rk_dt1.total_energy_density()(*s.m),
      s.rk_dt1.momentum_energy_density()(*s.m),
      s.rad.dt_radiation_energy_density_1(*s.m));

    actions::advance_rk_stage_1(s,
      sc,
      std::vector{std::make_tuple(s.rad.cons.radiation_energy_density(*s.m),
        s.rad.dt_radiation_energy_density_1(*s.m),
        s.rad.f.erad_face(s.m))},
      [&](auto & axis) {
        sc.execute<tasks::rad::compute_interface_fluxes<D>>(flecsi::exec::on,
          axis,
          *s.m,
          s.f.hydro.u_face(s.m),
          s.f.hydro.c_face(s.m),
          s.rad.f.erad_face(s.m),
          // Riemann Fluxes
          s.rad.rf.erad_f(*s.m),

          s.rad.dt_radiation_energy_density_1(*s.m));
      });
  }
};

static const auto rk_advance_1_action =
  spec::register_action<control, state, hydro_rk_advance_1, cp::rk_stage_1>();

template<std::size_t D>
struct update_1_vars {
  static void action(control_policy<state, D> & cp) {

    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::update_rk_stage_1(s,
      sc,
      std::vector{std::make_tuple(s.rad.cons.radiation_energy_density(*s.m),
        s.rad.dt_radiation_energy_density_1(*s.m),
        s.rad.dt_radiation_energy_density_n(*s.m))});
  }
};

static auto rk_vars_1_action =
  spec::register_action<control, state, update_1_vars, cp::rk_stage_1>();
static const bool rk_1_dp =
  spec::add_dependency(rk_vars_1_action, rk_advance_1_action);

} // namespace hard
