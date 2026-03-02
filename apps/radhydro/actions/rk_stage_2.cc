#include "rad/tasks/utils.hh"
#include "state.hh"

#include "../../actions/rk_stages/rk_stage_2.hh"

#include "rad/tasks/interface_fluxes.hh"
#include "rad/tasks/rad.hh"

namespace hard {

template<std::size_t D>
void
RK_advance_2(control_policy<state, D> & cp) {

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
      s.rad.f.erad_face(s.m))});

  // Calculate K1 and save it to dt_U
  // sc.execute<tasks::rad::compute_interface_fluxes<D>>(flecsi::exec::on,
  //  axis,
  //  *s.m,
  //  s.f.hydro.u_face(s.m),
  //  s.f.hydro.c_face(s.m),
  //  s.rad.f.erad_face(s.m),
  // Riemann Fluxes
  //  s.rad.rf.erad_f(*s.m),
  //
  //  s.rad.dt_radiation_energy_density_1(*s.m));
  //}
}

template<std::size_t D>
void
update_vars_2(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  actions::update_rk_stage_2(s,
    sc,
    std::vector{std::make_tuple(s.rad.cons.radiation_energy_density(*s.m),
      s.rad.dt_radiation_energy_density_1(*s.m),
      s.rad.dt_radiation_energy_density_2(*s.m),
      s.rad.dt_radiation_energy_density_n(*s.m))});

} // update_vars

inline control<state, 1>::action<RK_advance_2<1>, cp::rk_stage_2> rk_stage_2_1d;
inline control<state, 2>::action<RK_advance_2<2>, cp::rk_stage_2> rk_stage_2_2d;
inline control<state, 3>::action<RK_advance_2<3>, cp::rk_stage_2> rk_stage_2_3d;

inline control<state, 1>::action<update_vars_2<1>, cp::rk_stage_2>
  rk_stage_2_update_1d;
inline control<state, 2>::action<update_vars_2<2>, cp::rk_stage_2>
  rk_stage_2_update_2d;
inline control<state, 3>::action<update_vars_2<3>, cp::rk_stage_2>
  rk_stage_2_update_3d;

inline const auto dep_update_2_1d = rk_stage_2_update_1d.add(rk_stage_2_1d);
inline const auto dep_update_2_2d = rk_stage_2_update_2d.add(rk_stage_2_2d);
inline const auto dep_update_2_3d = rk_stage_2_update_3d.add(rk_stage_2_3d);

} // namespace hard
