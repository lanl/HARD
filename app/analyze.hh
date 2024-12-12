#ifndef FLASTRO_ANALYZE_HH
#define FLASTRO_ANALYZE_HH

#include "tasks/io.hh"
#include "types.hh"

#include <flecsi/flog.hh>
#include <spec/io.hh>

#ifdef USE_CATALYST
#include "tasks/catalyst.hh"
#endif

namespace flastro::action {

template<std::size_t D>
void
analyze(control_policy<state, D> & cp) {
  using namespace flecsi;
  auto & s = cp.state();
  auto lm = data::launch::make(s.m);

#ifndef FLASTRO_BENCHMARK_MODE
  if(((cp.step() % cp.output_frequency()) == 0) or
     (cp.step() == cp.max_steps())) {

    execute<tasks::io::raw<D>, mpi>(spec::io::name{"output-"}
                                      << std::setfill('0') << std::setw(5)
                                      << cp.step(),
      s.t(s.gt),
      lm,
      s.mass_density(lm),
      s.pressure(lm),
      s.velocity(lm),
      s.momentum_density(lm),
      s.total_energy_density(lm),
      s.radiation_energy_density(lm));

#ifdef USE_CATALYST

    if constexpr(D == 3) {
      // First update fields (catalyst_attributes) for catalysts pipeline
      flog(info) << "analyze action: prepare catalyst data" << std::endl;
      execute<tasks::external::update_attributes<D>, mpi>(catalyst_data(pt),
        s.t(s.gt),
        s.m,
        s.mass_density(s.m),
        s.velocity(s.m),
        s.pressure(s.m),
        s.total_energy_density(s.m),
        s.radiation_energy_density(s.m)); // <<< add variables here for catalyst

      // Then pass catalyst_data to catalyst pipeline
      flog(info) << "analyze action: execute catalyst" << std::endl;
      execute<tasks::external::execute_catalyst, mpi>(
        catalyst_data(pt), cp.step(), s.t(s.gt), lattice);
    }
    else {
      /* Do nothing, catalyst/paraview visualization only for 3D */
    }

#endif
  } // if
#endif

} // analyze

inline control<state, 1>::action<analyze<1>, cp::analyze> analyze1_action;
inline control<state, 2>::action<analyze<2>, cp::analyze> analyze2_action;
inline control<state, 3>::action<analyze<3>, cp::analyze> analyze3_action;

} // namespace flastro::action

#endif // FLASTRO_ANALYZE_HH
