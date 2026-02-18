#ifndef HARD_HYDRO_TIME_STEP_SIZE_HH
#define HARD_HYDRO_TIME_STEP_SIZE_HH

#include "state.hh"

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
time_step_size(control_policy<state, D> & cp) {
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

inline control<state, 1>::action<time_step_size<1>, cp::update_time_step_size>
  utss_1d;
inline control<state, 2>::action<time_step_size<2>, cp::update_time_step_size>
  utss_2d;
inline control<state, 3>::action<time_step_size<3>, cp::update_time_step_size>
  utss_3d;

#endif // HARD_HYDRO_TIME_STEP_SIZE_HH
