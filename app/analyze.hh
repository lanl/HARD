#ifndef HARD_ANALYZE_HH
#define HARD_ANALYZE_HH

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
analyze(control_policy<state, D> & cp) {
  using namespace flecsi;
  auto & s = cp.state();
  auto & sc = cp.scheduler();
  auto lm = data::launch::make(sc, *s.m);

#ifndef HARD_BENCHMARK_MODE

#if FLECSI_BACKEND == FLECSI_BACKEND_legion
  if(((cp.step() % cp.output_frequency()) == 0) or
     (cp.step() == cp.max_steps())) {
#else
  if(((cp.step() % cp.output_frequency()) == 0) or
     (cp.step() == cp.max_steps()) or (cp.time() == cp.max_time())) {
#endif
    execute<tasks::io::raw<D>, mpi>(flecsi::exec::on,
      spec::io::name{"output-"} << std::setfill('0') << std::setw(5)
                                << cp.step(),
      s.t(*s.gt),
      lm,
      s.mass_density(lm),
      s.pressure(lm),
      s.sound_speed(lm),
      s.specific_internal_energy(lm),
      s.velocity(lm),
      s.momentum_density(lm),
      s.total_energy_density(lm),
      s.radiation_energy_density(lm));

#ifdef USE_CATALYST

    if constexpr(D == 3) {
      // First update fields (catalyst_attributes) for catalysts pipeline
      flog(info) << "analyze action: prepare catalyst data" << std::endl;
      execute<tasks::external::update_attributes<D>, mpi>(flecsi::exec::on,
        s.catalyst_data(*s.pt),
        s.t(*s.gt),
        *s.m,
        s.mass_density(*s.m),
        s.velocity(*s.m),
        s.pressure(*s.m),
        s.total_energy_density(*s.m),
        s.radiation_energy_density(
          *s.m)); // <<< add variables here for catalyst

      // Then pass catalyst_data to catalyst pipeline
      flog(info) << "analyze action: execute catalyst" << std::endl;
      execute<tasks::external::execute_catalyst, mpi>(flecsi::exec::on,
        s.catalyst_data(*s.pt),
        cp.step(),
        s.t(*s.gt),
        lattice);
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

} // namespace hard::action

#endif // HARD_ANALYZE_HH
