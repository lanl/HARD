#include "state.hh"

#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/interface_fluxes.hh"
#include "hydro/tasks/reconstruct.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard {

template<std::size_t D>
void
update_ls(control_policy<state, D> & cp) {

  std::cout << "Hello World" << std::endl;

} // update_vars

inline control<state, 1>::action<update_ls<1>, cp::update_ls> update_ls_1d;
inline control<state, 2>::action<update_ls<2>, cp::update_ls> update_ls_2d;
inline control<state, 3>::action<update_ls<3>, cp::update_ls> update_ls_3d;

} // namespace hard
