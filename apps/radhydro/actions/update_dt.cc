#include "rad/tasks/utils.hh"
#include "state.hh"

#include "../../actions/update_dt/update_dt.hh"

namespace hard {

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
update_dt(control_policy<state, D> & cp) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  actions::update_dt(s, sc);

#ifdef HARD_ENABLE_LEGION_TRACING
  cp.guard.reset();
#endif
} // update_time_step_size

inline control<state, 1>::action<update_dt<1>, cp::update_dt> udt_1d;
inline control<state, 2>::action<update_dt<2>, cp::update_dt> udt_2d;
inline control<state, 3>::action<update_dt<3>, cp::update_dt> udt_3d;

} // namespace hard
