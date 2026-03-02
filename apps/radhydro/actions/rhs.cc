#include "rad/tasks/utils.hh"
#include "state.hh"

#include "../../actions/rhs/rhs.hh"

namespace hard {

template<std::size_t D>
void
rhs(control_policy<state, D> & cp) {
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

inline control<state, 1>::action<rhs<1>, cp::rhs> rhs_1d;
inline control<state, 2>::action<rhs<2>, cp::rhs> rhs_2d;
inline control<state, 3>::action<rhs<3>, cp::rhs> rhs_3d;

} // namespace hard
