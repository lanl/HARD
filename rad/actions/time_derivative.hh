#ifndef HARD_HYDRO_TIME_DERIVATIVE_HH
#define HARD_HYDRO_TIME_DERIVATIVE_HH

#include "state.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/init.hh"
#include "../modules/hydro/tasks/time_derivative.hh"

namespace hard {

template<std::size_t D>
void time_derivative(control_policy<state, D> & cp);

} // namespace hard

#endif // HARD_HYDRO_TIME_DERIVATIVE_HH
