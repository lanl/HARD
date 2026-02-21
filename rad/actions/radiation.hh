#ifndef HARD_RAD_RADIATION_HH
#define HARD_RAD_RADIATION_HH

#include "state.hh"

#include "../modules/hydro/tasks/cons2prim.hh"
#include "../modules/rad/linsolve.hh"
#include "../modules/rad/tasks/rad_root.hh"

namespace hard {

template<std::size_t D>

void radiation(control_policy<state, D> &);

inline control<state, 1>::action<radiation<1>, cp::radiation> radiation_1d;
inline control<state, 2>::action<radiation<2>, cp::radiation> radiation_2d;
inline control<state, 3>::action<radiation<3>, cp::radiation> radiation_3d;

} // namespace hard

#endif // HARD_RAD_RADIATION_HH
