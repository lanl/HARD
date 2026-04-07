#include "state.hh"

#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/interface_fluxes.hh"
#include "hydro/tasks/reconstruct.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

#include "../../actions/rk_stages/rk_stage_1.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct hydro_rk_advance_1 {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::advance_rk_stage_1(s, sc);
  }
};

static const auto rk_advance_1_action =
  spec::register_action<control, state, hydro_rk_advance_1, cp::rk_stage_1>();

template<std::size_t D>
struct hydro_update_1_vars {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::update_rk_stage_1(s, sc);
  }
};

static auto rk_vars_1_action =
  spec::register_action<control, state, hydro_update_1_vars, cp::rk_stage_1>();
static const bool rk_1_dp =
  spec::add_dependency(rk_vars_1_action, rk_advance_1_action);

} // namespace hard
