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

  actions::analyze<D>(cp,
    std::vector<std::tuple<actions::field_multi_double<D>, std::string>>{
      std::make_tuple(s.rad.cons.radiation_energy_density(lm),
        std::string{"radiation_energy_density"})});
} // analyze

inline control<state, 1>::action<analyze<1>, cp::initial_analyze>
  initial_analyze_1d;
inline control<state, 2>::action<analyze<2>, cp::initial_analyze>
  initial_analyze_2d;
inline control<state, 3>::action<analyze<3>, cp::initial_analyze>
  initial_analyze_3d;

inline control<state, 1>::action<analyze<1>, cp::analyze> analyze_1d;
inline control<state, 2>::action<analyze<2>, cp::analyze> analyze_2d;
inline control<state, 3>::action<analyze<3>, cp::analyze> analyze_3d;

} // namespace hard
