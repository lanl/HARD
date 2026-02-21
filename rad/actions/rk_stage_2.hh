#ifndef HARD_RAD_RK_STAGE_2_HH
#define HARD_RAD_RK_STAGE_2_HH

#include "state.hh"
#include "utils.hh"

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
void RK_advance_2(control_policy<state, D> & cp);

template<std::size_t D>
void update_vars_2(control_policy<state, D> & cp);

} // namespace hard

#endif // HARD_RAD_RK_STAGE_HH
