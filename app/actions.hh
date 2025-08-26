
#pragma once

#include "app/types.hh"
#include "numerical_algorithms/time_stepper.hh"
#include "rad.hh"
#include "state.hh"
#include "tasks/boundary.hh"
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

namespace hard::actions {

template<std::size_t D>
void
initialize_time_derivative(control_policy<state, D> & cp) {
  auto & s = cp.state();

#ifdef HARD_ENABLE_LEGION_TRACING
  // Legion tracing: Skip first iteration
  if(cp.step() == 0)
    cp.tracing.skip();
  // Legion tracing: Create new guard
  cp.guard.emplace(cp.tracing);
#endif

  flecsi::execute<tasks::init::compute_dt_weighted>(
    s.dt(s.gt), s.dt_weighted(s.gt), hard::time_stepper::time_stepper_gamma);

  // Set all dU_dt temporaries to zero before adding time derivative terms
  flecsi::execute<tasks::set_dudt_to_zero<D>, flecsi::default_accelerator>(s.m,
    s.dt_mass_density(s.m),
    s.dt_momentum_density(s.m),
    s.dt_total_energy_density(s.m),
    s.dt_radiation_energy_density(s.m));
  flecsi::execute<tasks::set_dudt_to_zero<D>, flecsi::default_accelerator>(s.m,
    s.dt_mass_density_2(s.m),
    s.dt_momentum_density_2(s.m),
    s.dt_total_energy_density_2(s.m),
    s.dt_radiation_energy_density_2(s.m));

  // Store the current state of evolved variables (U^n) before performing a time
  // step
  flecsi::execute<tasks::store_current_state<D>, flecsi::default_accelerator>(
    s.m,
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.radiation_energy_density(s.m),
    //
    s.mass_density_n(s.m),
    s.momentum_density_n(s.m),
    s.total_energy_density_n(s.m),
    s.radiation_energy_density_n(s.m));
}

// --------------------------------------------------------------------
//              EXPLICIT PARTS OF EVOLUTION
// --------------------------------------------------------------------

template<std::size_t D>
void
RK_advance(control_policy<state, D> & cp, time_stepper::rk_stage Stage) {

  auto & s = cp.state();

#ifdef ENABLE_RADIATION
  flecsi::execute<task::rad::getGradE<D>, flecsi::default_accelerator>(
    s.m, s.radiation_energy_density(s.m), s.gradient_rad_energy(s.m));

  flecsi::execute<task::rad::getLambda<D>, flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.radiation_energy_density(s.m),
    s.gradient_rad_energy(s.m),
    s.magnitude_gradient_rad_energy(s.m),
    s.R_value(s.m),
    s.lambda_bridge(s.m),
    kappa(s.gt));

  flecsi::execute<task::rad::getTensorP<D>, flecsi::default_accelerator>(s.m,
    s.radiation_pressure_tensor(s.m),
    s.radiation_energy_density(s.m),
    s.gradient_rad_energy(s.m),
    s.magnitude_gradient_rad_energy(s.m),
    s.lambda_bridge(s.m),
    s.R_value(s.m));

  flecsi::execute<task::rad::getRadForce<D>, flecsi::default_accelerator>(s.m,
    s.lambda_bridge(s.m),
    s.gradient_rad_energy(s.m),
    s.radiation_force(s.m));

  flecsi::execute<task::rad::getGradV<D>, flecsi::default_accelerator>(
    s.m, s.velocity_gradient(s.m), s.velocity(s.m));

  if(Stage == time_stepper::rk_stage::First) {
    flecsi::execute<task::rad::explicitSourceUpdate<D>,
      flecsi::default_accelerator>(s.m,
      s.velocity(s.m),
      s.radiation_force(s.m),
      s.radiation_pressure_tensor(s.m),
      s.velocity_gradient(s.m),
      //
      s.dt_momentum_density(s.m),
      s.dt_total_energy_density(s.m),
      s.dt_radiation_energy_density(s.m));
  }
  else if(Stage == time_stepper::rk_stage::Second) {
    flecsi::execute<task::rad::explicitSourceUpdate<D>,
      flecsi::default_accelerator>(s.m,
      s.velocity(s.m),
      s.radiation_force(s.m),
      s.radiation_pressure_tensor(s.m),
      s.velocity_gradient(s.m),
      //
      s.dt_momentum_density_2(s.m),
      s.dt_total_energy_density_2(s.m),
      s.dt_radiation_energy_density_2(s.m));
  }
#endif

  //
  // Perform reconstruction on cell faces, compute face fluxes with a Riemann
  // solver, and add the summation of dF^i/dx^i into (dU_dt)_explicit
  //

  using limiter = spec::limiters::weno5z;

  for(std::size_t axis = 0; axis < D; axis++) {
    // clang-format off
    flecsi::execute<tasks::hydro::reconstruct<D, limiter>,
      flecsi::default_accelerator>(axis,
      s.m, s.mass_density(s.m), s.velocity(s.m), s.pressure(s.m),
      s.specific_internal_energy(s.m),
      s.sound_speed(s.m),
      s.radiation_energy_density(s.m),
      s.rTail(s.m), s.rHead(s.m), s.uTail(s.m), s.uHead(s.m),
      s.pTail(s.m), s.pHead(s.m),  s.eTail(s.m), s.eHead(s.m),  s.cTail(s.m),
      s.cHead(s.m), s.EradTail(s.m), s.EradHead(s.m),
      s.ruTail(s.m), s.ruHead(s.m), s.rETail(s.m), s.rEHead(s.m));

    if(Stage == time_stepper::rk_stage::First) {
      // Calculate K1 and save it to dt_U
      flecsi::execute<tasks::hydro::compute_interface_fluxes<D>,
        flecsi::default_accelerator>(axis, s.m,
        s.rTail(s.m), s.rHead(s.m), s.uTail(s.m), s.uHead(s.m),
        s.pTail(s.m), s.pHead(s.m), s.eTail(s.m), s.eHead(s.m), s.cTail(s.m),
        s.cHead(s.m), s.EradTail(s.m), s.EradHead(s.m),
        s.ruTail(s.m), s.ruHead(s.m),
        s.rETail(s.m), s.rEHead(s.m),
        s.rF(s.m), s.ruF(s.m), s.rEF(s.m), s.EradF(s.m),
        s.dt_mass_density(s.m),
        s.dt_momentum_density(s.m),
        s.dt_total_energy_density(s.m),
        s.dt_radiation_energy_density(s.m));
    }
    else if(Stage == time_stepper::rk_stage::Second) {
      // Calculate K2 and save it to dt_U_2
      flecsi::execute<tasks::hydro::compute_interface_fluxes<D>,
        flecsi::default_accelerator>(axis, s.m,
        s.rTail(s.m), s.rHead(s.m), s.uTail(s.m), s.uHead(s.m),
        s.pTail(s.m), s.pHead(s.m), s.eTail(s.m), s.eHead(s.m), s.cTail(s.m),
        s.cHead(s.m), s.EradTail(s.m), s.EradHead(s.m),
        s.ruTail(s.m), s.ruHead(s.m),
        s.rETail(s.m), s.rEHead(s.m),
        s.rF(s.m), s.ruF(s.m), s.rEF(s.m), s.EradF(s.m),
        s.dt_mass_density_2(s.m),
        s.dt_momentum_density_2(s.m),
        s.dt_total_energy_density_2(s.m),
        s.dt_radiation_energy_density_2(s.m));
    }
    // clang-format on
  }
} // RK_advance

template<std::size_t D>
void
update_vars(control_policy<state, D> & cp, time_stepper::rk_stage Stage) {

  auto & s = cp.state();

  if(Stage == time_stepper::rk_stage::Second) {
    // Apply K1 to U with a Forward Euler step, so we can use U for the
    // K2 calculation in the next RK advance
    flecsi::execute<tasks::update_u<D>, flecsi::default_accelerator>(s.dt(s.gt),
      s.m,
      //
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      //
      s.dt_mass_density(s.m),
      s.dt_momentum_density(s.m),
      s.dt_total_energy_density(s.m),
      s.dt_radiation_energy_density(s.m));
  }
  else if(Stage == time_stepper::rk_stage::Update) {
    // First compute K1' = (K1 + K2) * 0.5
    flecsi::execute<tasks::add_k1_k2<D>, flecsi::default_accelerator>(s.m,
      //
      s.dt_mass_density(s.m),
      s.dt_momentum_density(s.m),
      s.dt_total_energy_density(s.m),
      s.dt_radiation_energy_density(s.m),
      //
      s.dt_mass_density_2(s.m),
      s.dt_momentum_density_2(s.m),
      s.dt_total_energy_density_2(s.m),
      s.dt_radiation_energy_density_2(s.m));

    // Now get U_n(+1) = U_n + h * K1'
    flecsi::execute<tasks::update_u<D>, flecsi::default_accelerator>(s.dt(s.gt),
      s.m,
      //
      s.mass_density_n(s.m),
      s.momentum_density_n(s.m),
      s.total_energy_density_n(s.m),
      s.radiation_energy_density_n(s.m),
      //
      s.dt_mass_density(s.m),
      s.dt_momentum_density(s.m),
      s.dt_total_energy_density(s.m),
      s.dt_radiation_energy_density(s.m));

    // Finish by updating the values stored in U_n to U
    flecsi::execute<tasks::store_current_state<D>, flecsi::default_accelerator>(
      s.m,
      s.mass_density_n(s.m),
      s.momentum_density_n(s.m),
      s.total_energy_density_n(s.m),
      s.radiation_energy_density_n(s.m),
      //
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m));
  }

  // Perform primitive recovery
  flecsi::execute<tasks::hydro::conservative_to_primitive<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.specific_internal_energy(s.m),
    s.sound_speed(s.m),
    s.eos);

  // Update boundary cells
  flecsi::execute<tasks::apply_boundaries<D>, flecsi::default_accelerator>(s.m,
    s.bmap(s.gt),
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.specific_internal_energy(s.m),
    s.radiation_energy_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m));

  if(s.mg) {
    // FIXME: figure out how not to use the hardcoded radiation temperature
    // boundary
    auto radiation_boundary_f =
      flecsi::execute<task::rad::interp_e_boundary>(s.t(s.gt),
        time_boundary(s.dense_topology),
        temperature_boundary(s.dense_topology));
    flecsi::execute<tasks::apply_radiation_boundary<D>,
      flecsi::default_accelerator>(
      s.m, s.radiation_energy_density(s.m), radiation_boundary_f);
  }

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

#ifdef ENABLE_RADIATION

  using namespace flecsi;
  auto & s = cp.state();

  flecsi::execute<task::rad_root::update_energy_density<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.velocity(s.m),
    s.temperature(s.m),
    s.total_energy_density(s.m),
    s.radiation_energy_density(s.m),
    kappa(s.gt),
    particle_mass(s.gt),
    s.dt_weighted(s.gt),
    s.eos);

  flecsi::execute<task::rad::getGradE<D>, flecsi::default_accelerator>(
    s.m, s.radiation_energy_density(s.m), s.gradient_rad_energy(s.m));

  flecsi::execute<task::rad::getLambda<D>, flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.radiation_energy_density(s.m),
    s.gradient_rad_energy(s.m),
    s.magnitude_gradient_rad_energy(s.m),
    s.R_value(s.m),
    s.lambda_bridge(s.m),
    kappa(s.gt));

  flecsi::execute<tasks::apply_boundary_single_field<D>,
    flecsi::default_accelerator>(s.m, s.bmap(s.gt), s.lambda_bridge(s.m));

  if(s.mg) {
    // FIXME: figure out how not to use the hardcoded radiation temperature
    // boundary
    auto radiation_boundary_f =
      flecsi::execute<task::rad::interp_e_boundary>(s.t(s.gt),
        time_boundary(s.dense_topology),
        temperature_boundary(s.dense_topology));
    flecsi::execute<tasks::apply_radiation_boundary<D>,
      flecsi::default_accelerator>(
      s.m, s.radiation_energy_density(s.m), radiation_boundary_f);
  }

  execute<task::rad::getDiff<D>, flecsi::default_accelerator>(
    s.m, s.mass_density(s.m), s.lambda_bridge(s.m), s.Diff(s.m), kappa(s.gt));

  // Initialize the diffusion coefficient
  execute<task::rad::diffusion_init<D>, flecsi::default_accelerator>(
    s.m, s.Diff(s.m), s.Df_x(s.m), s.Df_y(s.m), s.Df_z(s.m));

  // Initialize the stencil
  execute<task::rad::stencil_init<D>, flecsi::default_accelerator>(
    s.m, s.Df_x(s.m), s.Df_y(s.m), s.Df_z(s.m), s.Ew(s.m), s.dt_weighted(s.gt));

  // Initialize fields
  flecsi::execute<task::rad::initialize_Ef<D>, flecsi::default_accelerator>(
    s.m, s.radiation_energy_density(s.m), s.Ef(s.m));
  flecsi::execute<task::rad::const_init<D>, flecsi::default_accelerator>(
    s.m, s.Esf(s.m), 0.0);
  flecsi::execute<task::rad::const_init<D>, flecsi::default_accelerator>(
    s.m, s.Esf(s.m, 1), 0.0);
  flecsi::execute<task::rad::const_init<D>, flecsi::default_accelerator>(
    s.m, s.Resf(s.m), 0.0);

  hard::rad::fmg<D>(s);

  // Perform primitive recovery, since energy densities have changed
  flecsi::execute<tasks::hydro::conservative_to_primitive<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.specific_internal_energy(s.m),
    s.sound_speed(s.m),
    s.eos);

  // and also update boundary cells
  flecsi::execute<tasks::apply_boundaries<D>, flecsi::default_accelerator>(s.m,
    s.bmap(s.gt),
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.specific_internal_energy(s.m),
    s.radiation_energy_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m));

#endif

} // radiation_advance

// -----------------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next
//  time step.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
update_time_step_size(control_policy<state, D> & cp) {
  auto & s = cp.state();

  auto lmax_f =
    flecsi::execute<tasks::hydro::update_max_characteristic_speed<D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.velocity(s.m),
      s.pressure(s.m),
      s.sound_speed(s.m));

  s.dtmin_ =
    flecsi::reduce<hard::task::rad::update_dtmin<D>, flecsi::exec::fold::min>(
      s.m, lmax_f);

#ifdef HARD_ENABLE_LEGION_TRACING
  cp.guard.reset();
#endif
} // update_time_step_size

} // namespace hard::actions
