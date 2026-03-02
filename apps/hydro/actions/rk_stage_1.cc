#include "state.hh"

#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/interface_fluxes.hh"
#include "hydro/tasks/reconstruct.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

#include "../../actions/rk_stages/rk_stage_1.hh"

namespace hard {

template<std::size_t D>
void
hydro_RK_advance_1(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  actions::advance_rk_stage_1(s, sc);
}

template<std::size_t D>
void
hydro_update_vars(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  actions::update_rk_stage_1(s, sc);

} // update_vars

inline control<state, 1>::action<hydro_RK_advance_1<1>, cp::rk_stage_1>
  hydro_rk_stage_1_1d;
inline control<state, 2>::action<hydro_RK_advance_1<2>, cp::rk_stage_1>
  hydro_rk_stage_1_2d;
inline control<state, 3>::action<hydro_RK_advance_1<3>, cp::rk_stage_1>
  hydro_rk_stage_1_3d;

inline control<state, 1>::action<hydro_update_vars<1>, cp::rk_stage_1>
  hydro_rk_stage_1_update_1d;
inline control<state, 2>::action<hydro_update_vars<2>, cp::rk_stage_1>
  hydro_rk_stage_1_update_2d;
inline control<state, 3>::action<hydro_update_vars<3>, cp::rk_stage_1>
  hydro_rk_stage_1_update_3d;

inline const auto hydro_dep_update_1d =
  hydro_rk_stage_1_update_1d.add(hydro_rk_stage_1_1d);
inline const auto hydro_dep_update_2d =
  hydro_rk_stage_1_update_2d.add(hydro_rk_stage_1_2d);
inline const auto hydro_dep_update_3d =
  hydro_rk_stage_1_update_3d.add(hydro_rk_stage_1_3d);

} // namespace hard
