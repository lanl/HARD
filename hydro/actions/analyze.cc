#include "state.hh"
#include "utils.hh"

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
      std::vector{s.cons.hydro.mass_density(lm),
        s.prim.pressure(lm),
        s.prim.sound_speed(lm),
        s.prim.specific_internal_energy(lm),
        s.cons.hydro.total_energy_density(lm)},
      std::vector{s.prim.velocity(lm), s.cons.hydro.momentum_density(lm)},
      std::vector<std::string>{"density",
        "pressure",
        "sound_speed",
        "specific_internal_energy",
        "total_energy_density"},
      std::vector<std::string>{"velocity", "momentum_density"});
  } // if
#endif

} // analyze

inline control<state, 1>::action<analyze<1>, cp::analyze> analyze_1d;
inline control<state, 2>::action<analyze<2>, cp::analyze> analyze_2d;
inline control<state, 3>::action<analyze<3>, cp::analyze> analyze_3d;

} // namespace hard
