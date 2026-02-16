#pragma once

#include "app/types.hh"
#include "numerical_algorithms/time_stepper.hh"
#include "rad.hh"
#include "state.hh"
#include "tasks/boundaries/boundary.hh"
#include "tasks/external_source.hh"
#include "tasks/hydro/compute_interface_fluxes.hh"
#include "tasks/hydro/cons2prim.hh"
#include "tasks/hydro/maxcharspeed.hh"
#include "tasks/hydro/reconstruct.hh"
#include "tasks/init.hh"
#include "tasks/rad.hh"
#include "tasks/rad_root.hh"
#include "tasks/time_derivative.hh"
#include <cstddef>
#include <spec/limiter.hh>

#include "linsolve.hh"

namespace hard::actions {

// --------------------------------------------------------------------
//              EXPLICIT PARTS OF EVOLUTION
// --------------------------------------------------------------------

template<std::size_t D>
void
RK_advance(control_policy<state, D> & cp, time_stepper::rk_stage Stage) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  //
  // Perform reconstruction on cell faces, compute face fluxes with a Riemann
  // solver, and add the summation of dF^i/dx^i into (dU_dt)_explicit
  //

} // RK_advance

template<std::size_t D>
void
update_vars(control_policy<state, D> & cp, time_stepper::rk_stage Stage) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  if(Stage == time_stepper::rk_stage::Second) {
    // Apply K1 to U with a Forward Euler step, so we can use U for the

    // K2 calculation in the next RK advance
    sc.execute<tasks::update_u_stage<D, time_stepper::rk_stage::First>>(
      flecsi::exec::on,
      s.dt(*s.gt),
      *s.m,
      s.rk_n(s.m),
      s.rk_dt1(s.m),
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.cons.rad.radiation_energy_density(*s.m));
  }
  else if(Stage == time_stepper::rk_stage::Update) {
    // First compute K1' = (K1 + K2) * 0.5
    sc.execute<tasks::add_k1_k2<D>>(
      flecsi::exec::on, *s.m, s.rk_dt1(s.m), s.rk_dt2(s.m));

    // Now get U_n(+1) = U_n + h * K1'
    sc.execute<tasks::update_u<D>>(
      flecsi::exec::on, s.dt(*s.gt), *s.m, s.rk_n(s.m), s.rk_dt1(s.m));

    // Finish by updating the values stored in U_n to U
    sc.execute<tasks::store_current_state<D>>(flecsi::exec::on,
      *s.m,
      s.rk_n(s.m),
      //
      std::tuple{s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.cons.rad.radiation_energy_density(*s.m),
        s.cons.hydro.momentum_density(*s.m)});
  }

  // Perform primitive recovery
  sc.execute<tasks::hydro::conservative_to_primitive<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.cons.hydro.momentum_density(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.pressure(*s.m),
    s.prim.specific_internal_energy(*s.m),
    s.prim.sound_speed(*s.m),
    s.eos);

  // Update boundary cells
  sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.cons.rad.radiation_energy_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m)},
    std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});

} // update_vars

template<std::size_t D>
void
advection_advance(control_policy<state, D> & cp) {

  // First RK advance (this gives us K1)
  RK_advance<D>(cp, time_stepper::rk_stage::First);

  // Update variable for second advance (this is u0 + h * K1)
  update_vars(cp, time_stepper::rk_stage::Second);

  // Second RK advance (this gives us K2)
  RK_advance<D>(cp, time_stepper::rk_stage::Second);

  // Final variable update
  update_vars(cp, time_stepper::rk_stage::Update);

} // advection_advance

// --------------------------------------------------------------------
//              IMPLICIT PARTS OF EVOLUTION
// --------------------------------------------------------------------

template<std::size_t D>
void
radiation_advance(control_policy<state, D> & cp) {

  using namespace flecsi;
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  sc.execute<task::rad_root::update_energy_density<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.temperature(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.cons.rad.radiation_energy_density(*s.m),
    s.icst.rad.kappa(*s.gt),
    s.dt_weighted(*s.gt),
    s.eos);

  sc.execute<task::rad::getGradE<D>>(flecsi::exec::on,
    *s.m,
    s.cons.rad.radiation_energy_density(*s.m),
    s.rad_limiter.gradient_rad_energy(*s.m));

  // Adaptive FLD Radiation Advance

  sc.execute<task::rad::getLambda<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.cons.rad.radiation_energy_density(*s.m),
    s.rad_limiter.gradient_rad_energy(*s.m),
    s.rad_limiter.magnitude_gradient_rad_energy(*s.m),
    s.rad_limiter.R_value(*s.m),
    s.rad_limiter.lambda_bridge(*s.m),
    s.icst.rad.kappa(*s.gt),
    s.icst.rad.limiter_id(*s.gt));

  sc.execute<tasks::apply_boundaries_scalar<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.rad_limiter.lambda_bridge(*s.m)});

  sc.execute<task::rad::getDiff<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.rad_limiter.lambda_bridge(*s.m),
    s.mgr.Diff(*s.m),
    s.icst.rad.kappa(*s.gt));

  // Initialize the diffusion coefficient
  sc.execute<task::rad::diffusion_init<D>>(flecsi::exec::on,
    *s.m,
    s.mgr.Diff(*s.m),
    s.mgr.Df_x(*s.m),
    s.mgr.Df_y(*s.m),
    s.mgr.Df_z(*s.m));

  // Initialize the stencil
  sc.execute<task::rad::stencil_init<D>>(flecsi::exec::on,
    *s.m,
    s.mgr.Df_x(*s.m),
    s.mgr.Df_y(*s.m),
    s.mgr.Df_z(*s.m),
    s.mgr.Ew(*s.m),
    s.dt_weighted(*s.gt));

  // Initialize fields
  sc.execute<task::rad::copy_field<D>>(flecsi::exec::on,
    *s.m,
    s.cons.rad.radiation_energy_density(*s.m),
    s.mgr.Ef(*s.m));
  sc.execute<task::rad::const_init<D>>(
    flecsi::exec::on, *s.m, s.mgr.Esf(*s.m, 1), 0.0);
  sc.execute<task::rad::const_init<D>>(
    flecsi::exec::on, *s.m, s.mgr.Resf(*s.m), 0.0);

  std::chrono::time_point<std::chrono::system_clock> start_timer_rad =
    std::chrono::system_clock::now();

  hard::rad::linsolve<D>(cp);

  std::chrono::time_point<std::chrono::system_clock> stop_timer_rad =
    std::chrono::system_clock::now();

  flog(info) << " Radiation Timing: "
             << (stop_timer_rad - start_timer_rad).count() * 1e-9 << " [s] "
             << std::endl;

  // Move solution from rad solver
  sc.execute<task::rad::copy_field<D>>(flecsi::exec::on,
    *s.m,
    s.mgr.Uf(*s.m),
    s.cons.rad.radiation_energy_density(*s.m));

  // Perform primitive recovery, since energy densities have changed
  sc.execute<tasks::hydro::conservative_to_primitive<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.cons.hydro.momentum_density(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.pressure(*s.m),
    s.prim.specific_internal_energy(*s.m),
    s.prim.sound_speed(*s.m),
    s.eos);

  // and also update boundary cells
  sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.cons.rad.radiation_energy_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m)},
    std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});
} // radiation_advance

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
update_time_step_size(control_policy<state, D> & cp) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  auto lmax_f = sc.execute<tasks::hydro::update_max_characteristic_speed<D>>(
    flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.sound_speed(*s.m));

  s.dtmin_ =
    sc.reduce<hard::task::rad::update_dtmin<D>, flecsi::exec::fold::min>(
      flecsi::exec::on, *s.m, lmax_f);

#ifdef HARD_ENABLE_LEGION_TRACING
  cp.guard.reset();
#endif
} // update_time_step_size

} // namespace hard::actions
