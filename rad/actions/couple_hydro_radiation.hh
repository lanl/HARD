#ifndef HARD_RAD_COUPLE_HYDRO_RADIATION_HH
#define HARD_RAD_COUPLE_HYDRO_RADIATION_HH

#include "state.hh"

#include "../modules/rad/tasks/rad.hh"

namespace hard {

template<std::size_t D>
void couple_hydro_radiation(control_policy<state, D> &);

inline control<state, 1>::action<couple_hydro_radiation<1>,
  cp::couple_hydro_radiation_1>
  couple_hydro_radiation_1_1d;
inline control<state, 2>::action<couple_hydro_radiation<2>,
  cp::couple_hydro_radiation_1>
  couple_hydro_radiation_1_2d;
inline control<state, 3>::action<couple_hydro_radiation<3>,
  cp::couple_hydro_radiation_1>
  couple_hydro_radiation_1_3d;

inline control<state, 1>::action<couple_hydro_radiation<1>,
  cp::couple_hydro_radiation_2>
  couple_hydro_radiation_2_1d;
inline control<state, 2>::action<couple_hydro_radiation<2>,
  cp::couple_hydro_radiation_2>
  couple_hydro_radiation_2_2d;
inline control<state, 3>::action<couple_hydro_radiation<3>,
  cp::couple_hydro_radiation_2>
  couple_hydro_radiation_2_3d;

} // namespace hard

#endif // HARD_RAD_COUPLE_HYDRO_RADIATION_HH
