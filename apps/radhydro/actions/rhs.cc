#include "rad/tasks/utils.hh"
#include "state.hh"

#include "../../actions/rhs/rhs.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct rhs {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

#ifdef HARD_ENABLE_LEGION_TRACING
    // Legion tracing: Skip first iteration
    if(cp.step() == 0)
      cp.tracing.skip();
    // Legion tracing: Create new guard
    cp.guard.emplace(cp.tracing);
#endif

    actions::rhs(s,
      sc,
      std::vector{std::make_tuple(s.rad.cons.radiation_energy_density(*s.m),
        s.rad.dt_radiation_energy_density_1(*s.m),
        s.rad.dt_radiation_energy_density_2(*s.m),
        s.rad.dt_radiation_energy_density_n(*s.m))});
  }
};

static const auto rhs_action =
  spec::register_action<control, state, rhs, cp::rhs>();

} // namespace hard
