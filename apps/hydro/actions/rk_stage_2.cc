#include "state.hh"

#include "../../actions/rk_stages/rk_stage_2.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct hydro_rk_advance_2 {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::advance_rk_stage_2(s, sc);
  }
};

static const auto rk_advance_2_action =
  spec::register_action<control, state, hydro_rk_advance_2, cp::rk_stage_2>();

template<std::size_t D>
struct hydro_update_vars_2 {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    actions::update_rk_stage_2(s, sc);
  }
};

static auto rk_vars_2_action =
  spec::register_action<control, state, hydro_update_vars_2, cp::rk_stage_2>();
static const bool rk_2_dp =
  spec::add_dependency(rk_vars_2_action, rk_advance_2_action);

} // namespace hard
