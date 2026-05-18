#include "state.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct finalize {
  static void action(control_policy<state, D> &) {} // finalize
};

static const auto finalize_action =
  spec::register_action<control, state, finalize, cp::finalize>();

} // namespace hard
