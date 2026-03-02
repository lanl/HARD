#include "state.hh"

#include <flecsi/flog.hh>
#include <spec/tasks/io.hh>

#include "../../actions/analyze/analyze.hh"

namespace hard {

template<std::size_t D>
void
analyze(control_policy<state, D> & cp) {
  actions::analyze(cp);
} // analyze

inline control<state, 1>::action<analyze<1>, cp::analyze> analyze_1d;
inline control<state, 2>::action<analyze<2>, cp::analyze> analyze_2d;
inline control<state, 3>::action<analyze<3>, cp::analyze> analyze_3d;

} // namespace hard
