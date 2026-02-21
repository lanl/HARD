#ifndef HARD_HYDRO_INIT_HH
#define HARD_HYDRO_INIT_HH

#include "options.hh"
#include "state.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/cons2prim.hh"
#include "../modules/hydro/tasks/init.hh"
#include "../modules/hydro/tasks/initial_data/all_initial_data.hh"
#include "../modules/hydro/tasks/maxcharspeed.hh"
#include "../modules/hydro/tasks/time_derivative.hh"
#include "../modules/spec/eos.hh"
#include "../modules/spec/tasks/boundaries/boundary.hh"
#include "../modules/spec/tasks/io.hh"

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

namespace hard {

template<std::size_t D>
void initialize(control_policy<state, D> & cp);

inline control<state, 1>::action<initialize<1>, cp::initialize> init_1d;
inline control<state, 2>::action<initialize<2>, cp::initialize> init_2d;
inline control<state, 3>::action<initialize<3>, cp::initialize> init_3d;

} // namespace hard

#endif // HARD_HYDRO_INIT_HH
