#include "state.hh"

#include <flecsi/flog.hh>
#include <spec/io.hh>

#include <spec/runtime.hh>

namespace hard::action {

template<std::size_t D>
struct finalize {
  static void action(control_policy<state, D> &) {} // finalize
};

static const auto finalize_action =
  spec::register_action<control, state, finalize, cp::finalize>();

} // namespace hard::action
