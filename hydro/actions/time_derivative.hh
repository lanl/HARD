#ifndef HARD_HYDRO_TIME_DERIVATIVE_HH
#define HARD_HYDRO_TIME_DERIVATIVE_HH

#include "state.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/init.hh"
#include "../modules/hydro/tasks/time_derivative.hh"

namespace hard {

template<std::size_t D>
void time_derivative(control_policy<state, D> & cp);

inline control<state, 1>::action<time_derivative<1>, cp::time_derivative>
  time_derivative_1d;
inline control<state, 2>::action<time_derivative<2>, cp::time_derivative>
  time_derivative_2d;
inline control<state, 3>::action<time_derivative<3>, cp::time_derivative>
  time_derivative_3d;

} // namespace hard

#endif // HARD_HYDRO_TIME_DERIVATIVE_HH
