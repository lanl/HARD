#ifndef HARD_HYDRO_FINALIZE_HH
#define HARD_HYDRO_FINALIZE_HH

#include "state.hh"

#include <../modules/spec/io.hh>
#include <flecsi/flog.hh>

namespace hard::action {

template<std::size_t D>
void
finalize([[maybe_unused]] control_policy<state, D> & cs) {
  using namespace flecsi;
} // finalize

inline control<state, 1>::action<finalize<1>, cp::finalize> finalize1_action;
inline control<state, 2>::action<finalize<2>, cp::finalize> finalize2_action;
inline control<state, 3>::action<finalize<3>, cp::finalize> finalize3_action;

} // namespace hard::action

#endif // HARD_HYDRO_FINALIZE_HH
