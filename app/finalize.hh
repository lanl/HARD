#ifndef HARD_FINALIZE_HH
#define HARD_FINALIZE_HH

#include "tasks/io.hh"
#include "types.hh"

#include <flecsi/flog.hh>
#include <spec/io.hh>

#ifdef USE_CATALYST
#include "tasks/catalyst.hh"
#endif

namespace hard::action {

template<std::size_t D>
void
finalize(control_policy<state, D> & cp) {
  using namespace flecsi;
  auto & s = cp.state();

#ifdef USE_CATALYST

  if constexpr(D == 3) {
    execute<tasks::external::finalize, mpi>(catalyst_data(pt));
    flog(info) << "finalize action: done" << std::endl;
  }
  else {
    /* Do nothing, catalyst/paraview visualization only for 3D */
  }

#endif

} // finalize

inline control<state, 1>::action<finalize<1>, cp::finalize> finalize1_action;
inline control<state, 2>::action<finalize<2>, cp::finalize> finalize2_action;
inline control<state, 3>::action<finalize<3>, cp::finalize> finalize3_action;

} // namespace hard::action

#endif // HARD_FINALIZE_HH
