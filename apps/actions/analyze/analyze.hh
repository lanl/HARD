#ifndef HARD_APPS_ACTIONS_ANALYZE_HH
#define HARD_APPS_ACTIONS_ANALYZE_HH

namespace hard::actions {

template<std::size_t D>
using field_multi_double = flecsi::data::
  multi_reference<double, flecsi::data::dense, spec::mesh<D>, spec::is::cells>;

template<std::size_t D>
void
analyze([[maybe_unused]] control_policy<state, D> & cp,
  [[maybe_unused]] std::vector<std::tuple<field_multi_double<D>, const char *>>
    v_t = {}) {

#ifndef HARD_BENCHMARK_MODE

#if FLECSI_BACKEND == FLECSI_BACKEND_legion
  if(((cp.step() % cp.output_frequency()) == 0) or
     (cp.step() == cp.max_steps())) {
#else
  if(((cp.step() % cp.output_frequency()) == 0) or
     (cp.step() == cp.max_steps()) or (cp.time() == cp.max_time())) {
#endif

    using namespace flecsi;
    auto & s = cp.state();
    auto & sc = cp.scheduler();
    auto lm = data::launch::make(sc, *s.m);

    auto scalar_v =
      std::vector{std::make_tuple(s.cons.hydro.mass_density(lm), "density"),
        std::make_tuple(s.prim.pressure(lm), "pressure"),
        std::make_tuple(s.prim.sound_speed(lm), "sound_speed"),
        std::make_tuple(
          s.prim.specific_internal_energy(lm), "specific_internal_energy"),
        std::make_tuple(
          s.cons.hydro.total_energy_density(lm), "total_energy_density")};

    scalar_v.insert(scalar_v.end(), v_t.begin(), v_t.end());

    if(s.output_method == spec::om::vti) {
      flecsi::execute<tasks::io::vtk<D>, flecsi::mpi>(flecsi::exec::on,
        spec::io::name{""} << std::setfill('0') << std::setw(5) << cp.step(),
        s.t(*s.gt),
        lm,
        scalar_v,
        std::vector{std::make_tuple(s.prim.velocity(lm), "velocity"),
          std::make_tuple(
            s.cons.hydro.momentum_density(lm), "momentum_density")});
    }
    else {
      flecsi::execute<tasks::io::csv<D>, flecsi::mpi>(flecsi::exec::on,
        spec::io::name{""} << std::setfill('0') << std::setw(5) << cp.step(),
        s.t(*s.gt),
        lm,
        scalar_v,
        std::vector{std::make_tuple(s.prim.velocity(lm), "velocity"),
          std::make_tuple(
            s.cons.hydro.momentum_density(lm), "momentum_density")});
    }
  } // if

#endif

} // analyze

} // namespace hard::actions

#endif // HARD_APPS_ACTIONS_ANALYZE_HH
