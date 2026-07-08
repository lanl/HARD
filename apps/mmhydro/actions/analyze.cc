#include "state.hh"

#include "gfm/tasks/utils.hh"

#include <flecsi/flog.hh>
#include <spec/tasks/io.hh>
#include <string>
#include <vector>

#include "../../actions/analyze/analyze.hh"

namespace hard {

template<std::size_t D>
void
analyze(control_policy<state, D> & cp) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  sc.execute<tasks::util::material_to_scalar_output<D>>(
    flecsi::exec::on, *s.m, s.ls.material(*s.m), s.ls.material_out(*s.m));

  auto lm = flecsi::data::launch::make(sc, *s.m);

  actions::analyze<D>(cp,
    std::vector<std::tuple<actions::field_multi_double<D>, std::string>>{
      std::make_tuple(s.ls.psi(lm), std::string{"psi"}),
      std::make_tuple(s.ls.material_out(lm), std::string{"material_id"})});
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
