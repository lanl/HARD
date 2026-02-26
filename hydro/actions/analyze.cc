#include "state.hh"

#include <../modules/spec/tasks/io.hh>
#include <flecsi/flog.hh>

namespace hard {

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

    execute<tasks::io::csv<D>, mpi>(flecsi::exec::on,
      spec::io::name{""} << std::setfill('0') << std::setw(5) << cp.step(),
      s.t(*s.gt),
      lm,
      std::vector{std::make_tuple(s.cons.hydro.mass_density(lm), "density"),
        std::make_tuple(s.prim.pressure(lm), "pressure"),
        std::make_tuple(s.prim.sound_speed(lm), "sound_speed"),
        std::make_tuple(
          s.prim.specific_internal_energy(lm), "specific_internal_energy"),
        std::make_tuple(
          s.cons.hydro.total_energy_density(lm), "total_energy_density")},
      std::vector{std::make_tuple(s.prim.velocity(lm), "velocity"),
        std::make_tuple(
          s.cons.hydro.momentum_density(lm), "momentum_density")});

  } // if
#endif

} // analyze

inline control<state, 1>::action<analyze<1>, cp::analyze> analyze_1d;
inline control<state, 2>::action<analyze<2>, cp::analyze> analyze_2d;
inline control<state, 3>::action<analyze<3>, cp::analyze> analyze_3d;

} // namespace hard
