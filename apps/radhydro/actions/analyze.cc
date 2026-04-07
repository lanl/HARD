#include "rad/tasks/utils.hh"
#include "state.hh"

#include <flecsi/flog.hh>
#include <spec/tasks/io.hh>

#include "../../actions/analyze/analyze.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct analyze {

  static void action(control_policy<state, D> & cp) {

    using namespace flecsi;
    auto & s = cp.state();
    auto & sc = cp.scheduler();
    auto lm = data::launch::make(sc, *s.m);

    actions::analyze(cp,
      std::vector{std::make_tuple(
        s.rad.cons.radiation_energy_density(lm), "radiation_energy_density")});
  }
};

static const auto analyze_action =
  spec::register_action<control, state, analyze, cp::analyze>();

} // namespace hard
