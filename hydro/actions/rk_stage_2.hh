#ifndef HARD_HYDRO_RK_STAGE_2_HH
#define HARD_HYDRO_RK_STAGE_2_HH

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
void RK_advance_2(control_policy<state, D> & cp);

template<std::size_t D>
void update_vars_2(control_policy<state, D> & cp);

inline control<state, 1>::action<RK_advance_2<1>, cp::rk_stage_2>
  hydro_rk_stage_2_1d;
inline control<state, 2>::action<RK_advance_2<2>, cp::rk_stage_2>
  hydro_rk_stage_2_2d;
inline control<state, 3>::action<RK_advance_2<3>, cp::rk_stage_2>
  hydro_rk_stage_2_3d;

inline control<state, 1>::action<update_vars_2<1>, cp::rk_stage_2>
  hydro_rk_stage_2_update_1d;
inline control<state, 2>::action<update_vars_2<2>, cp::rk_stage_2>
  hydro_rk_stage_2_update_2d;
inline control<state, 3>::action<update_vars_2<3>, cp::rk_stage_2>
  hydro_rk_stage_2_update_3d;

inline const auto dep_update_2_1d =
  hydro_rk_stage_2_update_1d.add(hydro_rk_stage_2_1d);
inline const auto dep_update_2_2d =
  hydro_rk_stage_2_update_2d.add(hydro_rk_stage_2_2d);
inline const auto dep_update_2_3d =
  hydro_rk_stage_2_update_3d.add(hydro_rk_stage_2_3d);

} // namespace hard

#endif // HARD_HYDRO_RK_STAGE_2_HH
