#ifndef HARD_RAD_INIT_HH
#define HARD_RAD_INIT_HH

#include "options.hh"
#include "state.hh"
#include "types.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/cons2prim.hh"
#include "../modules/hydro/tasks/init.hh"
#include "../modules/hydro/tasks/maxcharspeed.hh"
#include "../modules/hydro/tasks/time_derivative.hh"
#include "../modules/rad/tasks/init.hh"
#include "../modules/rad/tasks/initial_data/all_initial_data.hh"
#include "../modules/spec/eos.hh"
#include "../modules/spec/tasks/boundaries/boundary.hh"
#include "../modules/spec/tasks/io.hh"

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

namespace hard {

template<std::size_t D>
void initialize(control_policy<state, D> & cp);

} // namespace hard

#endif // HARD_HYDRO_INITIALIZE_HH
