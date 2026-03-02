#include "rad/tasks/utils.hh"
#include "state.hh"

#include <flecsi/flog.hh>
#include <spec/tasks/io.hh>

#include "../../actions/analyze/analyze.hh"

namespace hard {

template<std::size_t D>
void
analyze(control_policy<state, D> & cp) {

  using namespace flecsi;
  auto & s = cp.state();
  auto & sc = cp.scheduler();
  auto lm = data::launch::make(sc, *s.m);

  actions::analyze(cp,
    std::vector{std::make_tuple(
      s.rad.cons.radiation_energy_density(lm), "radiation_energy_density")});

} // analyze

inline control<state, 1>::action<analyze<1>, cp::analyze> analyze_1d;
inline control<state, 2>::action<analyze<2>, cp::analyze> analyze_2d;
inline control<state, 3>::action<analyze<3>, cp::analyze> analyze_3d;

} // namespace hard
