#include "state.hh"

#include "gfm/tasks/ls_core.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard {

template<std::size_t D>
void
update_ls(control_policy<state, D> & cp) {

  hard::tasks::level_set::ssprk3_step(cp);

} // update_ls

inline control<state, 1>::action<update_ls<1>, cp::update_ls> update_ls_1d;
inline control<state, 2>::action<update_ls<2>, cp::update_ls> update_ls_2d;
inline control<state, 3>::action<update_ls<3>, cp::update_ls> update_ls_3d;

} // namespace hard
