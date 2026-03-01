#include "rad/tasks/utils.hh"
#include "state.hh"

#include "hydro/tasks/init.hh"
#include "hydro/tasks/rhs.hh"

namespace hard {

template<std::size_t D>
void
rhs(control_policy<state, D> & cp) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

#ifdef HARD_ENABLE_LEGION_TRACING
  // Legion tracing: Skip first iteration
  if(cp.step() == 0)
    cp.tracing.skip();
  // Legion tracing: Create new guard
  cp.guard.emplace(cp.tracing);
#endif

  sc.execute<tasks::init::compute_dt_weighted>(flecsi::exec::on,
    s.dt(*s.gt),
    s.dt_weighted(*s.gt),
    hard::time_stepper::time_stepper_gamma);

  // Set all dU_dt temporaries to zero before adding time derivative terms
  sc.execute<tasks::hydro::set_dudt_to_zero<D>>(flecsi::exec::on,
    std::vector{// dt1
      s.rk_dt1.mass_density()(*s.m),
      s.rk_dt1.total_energy_density()(*s.m),
      s.rad.dt_radiation_energy_density_1(*s.m),
      // dt2
      s.rk_dt2.mass_density()(*s.m),
      s.rk_dt2.total_energy_density()(*s.m),
      s.rad.dt_radiation_energy_density_2(*s.m)},
    std::vector{// dt1
      s.rk_dt1.momentum_energy_density()(*s.m),
      // dt2
      s.rk_dt2.momentum_energy_density()(*s.m)});

  // Store the current state of evolved variables (U^n) before performing a time
  // step
  // clang-format off
  sc.execute<tasks::hydro::store_current_state<D>>(flecsi::exec::on,
    std::vector{
      std::make_tuple(
        s.cons.hydro.mass_density(*s.m),
        s.rk_n.mass_density()(*s.m)),
      std::make_tuple(
        s.cons.hydro.total_energy_density(*s.m),
        s.rk_n.total_energy_density()(*s.m)),
      std::make_tuple(
        s.rad.cons.radiation_energy_density(*s.m),
        s.rad.dt_radiation_energy_density_n(*s.m))},
    std::vector{
      std::make_tuple(
        s.cons.hydro.momentum_density(*s.m),
        s.rk_n.momentum_energy_density()(*s.m))});
  // clang-format on
}

inline control<state, 1>::action<rhs<1>, cp::rhs> rhs_1d;
inline control<state, 2>::action<rhs<2>, cp::rhs> rhs_2d;
inline control<state, 3>::action<rhs<3>, cp::rhs> rhs_3d;

} // namespace hard
