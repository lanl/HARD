#ifndef HARD_ADVANCE_HH
#define HARD_ADVANCE_HH

#include "actions.hh"
#include "state.hh"

#include <flecsi/execution.hh>

namespace hard::actions {

// ----------------------------------------------------------------------------
//              Action lists and dependencies
// ----------------------------------------------------------------------------

using namespace hard::time_stepper;

// Solve by operator split:
// 0) Initialize all variables that will be used
//
// 1) Advance radiation terms
//
// 2) Do radiation advection and hydro using Heun's method

inline control<state, 1>::action<initialize_time_derivative<1>, cp::advance>
  initialize_time_derivative_1d;
inline control<state, 2>::action<initialize_time_derivative<2>, cp::advance>
  initialize_time_derivative_2d;
inline control<state, 3>::action<initialize_time_derivative<3>, cp::advance>
  initialize_time_derivative_3d;

#ifdef ENABLE_RADIATION
// ----------------------------------------------------------------------------
//              Solve radiation diffusion and source-sink terms
// ----------------------------------------------------------------------------

inline control<state, 1>::action<radiation_advance<1>, cp::advance>
  radiation_advance_1d;
inline control<state, 2>::action<radiation_advance<2>, cp::advance>
  radiation_advance_2d;
inline control<state, 3>::action<radiation_advance<3>, cp::advance>
  radiation_advance_3d;
inline const auto dep_radiation_advance_1d =
  radiation_advance_1d.add(initialize_time_derivative_1d);
inline const auto dep_radiation_advance_2d =
  radiation_advance_2d.add(initialize_time_derivative_2d);
inline const auto dep_radiation_advance_3d =
  radiation_advance_3d.add(initialize_time_derivative_3d);
#endif

// ----------------------------------------------------------------------------
//              Solve advection
// ----------------------------------------------------------------------------

inline control<state, 1>::action<advection_advance<1>, cp::advance>
  advection_advance_1d;
inline control<state, 2>::action<advection_advance<2>, cp::advance>
  advection_advance_2d;
inline control<state, 3>::action<advection_advance<3>, cp::advance>
  advection_advance_3d;
#ifdef ENABLE_RADIATION
inline const auto dep_advection_advance_1d =
  advection_advance_1d.add(radiation_advance_1d);
inline const auto dep_advection_advance_2d =
  advection_advance_2d.add(radiation_advance_2d);
inline const auto dep_advection_advance_3d =
  advection_advance_3d.add(radiation_advance_3d);
#else
inline const auto dep_advection_advance_1d =
  advection_advance_1d.add(initialize_time_derivative_1d);
inline const auto dep_advection_advance_2d =
  advection_advance_2d.add(initialize_time_derivative_2d);
inline const auto dep_advection_advance_3d =
  advection_advance_3d.add(initialize_time_derivative_3d);
#endif
// --------------------------------------------------------------------
//              Update for next cycle
// --------------------------------------------------------------------

// Update the CFL limit
inline control<state, 1>::action<update_time_step_size<1>, cp::advance>
  update_time_step_size_1d;
inline control<state, 2>::action<update_time_step_size<2>, cp::advance>
  update_time_step_size_2d;
inline control<state, 3>::action<update_time_step_size<3>, cp::advance>
  update_time_step_size_3d;
inline const auto dep_update_time_step_size_1d =
  update_time_step_size_1d.add(advection_advance_1d);
inline const auto dep_update_time_step_size_2d =
  update_time_step_size_2d.add(advection_advance_2d);
inline const auto dep_update_time_step_size_3d =
  update_time_step_size_3d.add(advection_advance_3d);

} // namespace hard::actions

#endif // HARD_ADVANCE_HH
