#include "state.hh"

#include <spec/runtime.hh>
#include <spec/tasks/io.hh>

#include <flecsi/flog.hh>

#include "../../actions/analyze/analyze.hh"

namespace hard {

template<std::size_t D>
struct analyze {
  static void action(control_policy<state, D> & cp) {
    actions::analyze(cp);
  }
};

static const auto analyze_action =
  spec::register_action<control, state, analyze, cp::analyze>();

} // namespace hard
