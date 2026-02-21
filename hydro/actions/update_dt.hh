#ifndef HARD_HYDRO_UPDATE_DT_HH
#define HARD_HYDRO_UPDATE_DT_HH

#include "state.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/maxcharspeed.hh"
#include "../modules/hydro/tasks/time_derivative.hh"

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
namespace hard {

template<std::size_t D>
void update_dt(control_policy<state, D> & cp);

inline control<state, 1>::action<update_dt<1>, cp::update_dt> udt_1d;
inline control<state, 2>::action<update_dt<2>, cp::update_dt> udt_2d;
inline control<state, 3>::action<update_dt<3>, cp::update_dt> udt_3d;

} // namespace hard

#endif // HARD_HYDRO_TIME_STEP_SIZE_HH
