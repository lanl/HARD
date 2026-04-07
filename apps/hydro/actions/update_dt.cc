#include "state.hh"

#include "../../actions/update_dt/update_dt.hh"

#include <spec/runtime.hh>

namespace hard {

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
template<std::size_t D>
struct update_dt {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::update_dt(s, sc);

#ifdef HARD_ENABLE_LEGION_TRACING
    cp.guard.reset();
#endif
  }
};

static const auto update_dt_action =
  spec::register_action<control, state, update_dt, cp::update_dt>();

} // namespace hard
