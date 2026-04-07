#include "rad/tasks/utils.hh"
#include "state.hh"

#include "../../actions/rk_stages/rk_stage_2.hh"

#include "rad/tasks/interface_fluxes.hh"
#include "rad/tasks/rad.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct rk_advance_2 {
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
      s.rk_dt2.total_energy_density()(*s.m),
      s.rk_dt2.momentum_energy_density()(*s.m),
      s.rad.dt_radiation_energy_density_2(*s.m));

    actions::advance_rk_stage_2(s,
      sc,
      std::vector{std::make_tuple(s.rad.cons.radiation_energy_density(*s.m),
        s.rad.dt_radiation_energy_density_2(*s.m),
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

          s.rad.dt_radiation_energy_density_2(*s.m));
      });
  }
};

static const auto rk_advance_2_action =
  spec::register_action<control, state, rk_advance_2, cp::rk_stage_2>();

template<std::size_t D>
struct update_vars_2 {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::update_rk_stage_2(s,
      sc,
      std::vector{std::make_tuple(s.rad.cons.radiation_energy_density(*s.m),
        s.rad.dt_radiation_energy_density_1(*s.m),
        s.rad.dt_radiation_energy_density_2(*s.m),
        s.rad.dt_radiation_energy_density_n(*s.m))});
  }
};

static auto rk_vars_2_action =
  spec::register_action<control, state, update_vars_2, cp::rk_stage_2>();
static const bool rk_2_dp =
  spec::add_dependency(rk_vars_2_action, rk_advance_2_action);

} // namespace hard
