#ifndef HARD_HYDRO_RK_STAGE_1_HH
#define HARD_HYDRO_RK_STAGE_1_HH

#include "state.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/cons2prim.hh"
#include "../modules/hydro/tasks/external_source.hh"
#include "../modules/hydro/tasks/interface_fluxes.hh"
#include "../modules/hydro/tasks/reconstruct.hh"
#include "../modules/hydro/tasks/time_derivative.hh"
#include "../modules/spec/limiter.hh"
#include "../modules/spec/tasks/boundaries/boundary.hh"

namespace hard {

template<std::size_t D>
void hydro_RK_advance_1(control_policy<state, D> & cp);

template<std::size_t D>
void hydro_update_vars(control_policy<state, D> & cp);

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

#endif // HARD_HYDRO_RK_STAGE_1_HH
