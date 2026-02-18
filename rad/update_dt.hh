#ifndef HARD_HYDRO_UPDATE_DT_HH
#define HARD_HYDRO_UPDATE_DT_HH

#include "state.hh"

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
update_dt(control_policy<state, D> & cp) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  auto lmax_f = sc.execute<tasks::hydro::update_max_characteristic_speed<D>>(
    flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.sound_speed(*s.m));

  s.dtmin_ = sc.reduce<tasks::hydro::update_dtmin<D>, flecsi::exec::fold::min>(
    flecsi::exec::on, *s.m, lmax_f);

#ifdef HARD_ENABLE_LEGION_TRACING
  cp.guard.reset();
#endif
} // update_time_step_size

inline control<state, 1>::action<update_dt<1>, cp::update_dt> udt_1d;
inline control<state, 2>::action<update_dt<2>, cp::update_dt> udt_2d;
inline control<state, 3>::action<update_dt<3>, cp::update_dt> udt_3d;

#endif // HARD_HYDRO_TIME_STEP_SIZE_HH
