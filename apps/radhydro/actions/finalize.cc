#include "state.hh"

#include <spec/io.hh>
#include <flecsi/flog.hh>

namespace hard::action {

template<std::size_t D>
void
finalize(control_policy<state, D> &) {}

inline control<state, 1>::action<finalize<1>, cp::finalize> finalize1_action;
inline control<state, 2>::action<finalize<2>, cp::finalize> finalize2_action;
inline control<state, 3>::action<finalize<3>, cp::finalize> finalize3_action;

} // namespace hard::action
