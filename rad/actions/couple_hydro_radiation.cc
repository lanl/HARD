#ifndef HARD_RAD_COUPLE_HYDRO_RADIATION_HH
#define HARD_RAD_COUPLE_HYDRO_RADIATION_HH

#include "state.hh"

#include "../modules/rad/tasks/rad.hh"

using namespace hard;

template<std::size_t D>
void couple_hydro_radiation(control_policy<state, D> &);

#endif // HARD_RAD_COUPLE_HYDRO_RADIATION_HH
