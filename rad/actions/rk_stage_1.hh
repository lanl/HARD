#ifndef HARD_RAD_RK_STAGE_1_HH
#define HARD_RAD_RK_STAGE_1_HH

#include "state.hh"
#include "utils.hh"

#include "../hydro/actions/rk_stage_1.hh"

#include "../modules/hydro/tasks/cons2prim.hh"
#include "../modules/hydro/tasks/external_source.hh"
#include "../modules/hydro/tasks/interface_fluxes.hh"
#include "../modules/hydro/tasks/reconstruct.hh"
#include "../modules/hydro/tasks/time_derivative.hh"
#include "../modules/rad/tasks/interface_fluxes.hh"
#include "../modules/rad/tasks/rad.hh"
#include "../modules/spec/limiter.hh"
#include "../modules/spec/tasks/boundaries/boundary.hh"

namespace hard {

template<std::size_t D>
void RK_advance_1(control_policy<state, D> & cp);

template<std::size_t D>
void update_vars(control_policy<state, D> & cp);

inline control<state, 1>::action<RK_advance_1<1>, cp::rk_stage_1>
  rad_rk_stage_1_1d;
inline control<state, 2>::action<RK_advance_1<2>, cp::rk_stage_1>
  rad_rk_stage_1_2d;
inline control<state, 3>::action<RK_advance_1<3>, cp::rk_stage_1>
  rad_rk_stage_1_3d;

inline control<state, 1>::action<update_vars<1>, cp::rk_stage_1>
  rad_rk_stage_1_update_1d;
inline control<state, 2>::action<update_vars<2>, cp::rk_stage_1>
  rad_rk_stage_1_update_2d;
inline control<state, 3>::action<update_vars<3>, cp::rk_stage_1>
  rad_rk_stage_1_update_3d;

//  hydro_rk_stage_1 -> rad_rk_stage_1
inline const auto dep_hydro_1d = rad_rk_stage_1_1d.add(hydro_rk_stage_1_1d);
inline const auto dep_hydro_2d = rad_rk_stage_1_2d.add(hydro_rk_stage_1_2d);
inline const auto dep_hydro_3d = rad_rk_stage_1_3d.add(hydro_rk_stage_1_3d);

// rad_rk_stage_1 -> hydro_update_var
inline const auto dep_hydro_update_1d =
  hydro_rk_stage_1_update_1d.add(rad_rk_stage_1_1d);
inline const auto dep_hydro_update_2d =
  hydro_rk_stage_1_update_2d.add(rad_rk_stage_1_2d);
inline const auto dep_hydro_update_3d =
  hydro_rk_stage_1_update_3d.add(rad_rk_stage_1_3d);

// hydro_update_var -> rad_update_var
inline const auto dep_update_1d =
  rad_rk_stage_1_update_1d.add(hydro_rk_stage_1_update_1d);
inline const auto dep_update_2d =
  rad_rk_stage_1_update_2d.add(hydro_rk_stage_1_update_2d);
inline const auto dep_update_3d =
  rad_rk_stage_1_update_3d.add(hydro_rk_stage_1_update_3d);

} // namespace hard

#endif // HARD_RAD_RK_STAGE_1_HH
