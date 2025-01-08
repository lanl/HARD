
#pragma once

#include "numerical_algorithms/time_stepper.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/hydro.hh"
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
    s.dt_radiation_energy_density(s.m),
    s.dt_total_energy_density_implicit(s.m),
    s.dt_radiation_energy_density_implicit(s.m));
  flecsi::execute<tasks::set_dudt_to_zero<D>, flecsi::default_accelerator>(s.m,
    s.dt_mass_density_2(s.m),
    s.dt_momentum_density_2(s.m),
    s.dt_total_energy_density_2(s.m),
    s.dt_radiation_energy_density_2(s.m),
    s.dt_total_energy_density_implicit_2(s.m),
    s.dt_radiation_energy_density_implicit_2(s.m));

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

template<std::size_t D, time_stepper::rk_stage Stage>
void
explicit_source_terms(control_policy<state, D> & cp) {

#ifndef DISABLE_RADIATION

  auto & s = cp.state();

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

  // Compute explicit source terms and store them into (dU_dt)_explicit
  if constexpr(Stage == time_stepper::rk_stage::First) {
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
  else if constexpr(Stage == time_stepper::rk_stage::Second) {
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
}

//
// Perform reconstruction on cell faces, compute face fluxes with a Riemann
// solver, and add the summation of dF^i/dx^i into (dU_dt)_explicit
//
template<std::size_t D, time_stepper::rk_stage Stage>
void
fluxes_terms(control_policy<state, D> & cp) {
  auto & s = cp.state();

  using limiter = spec::limiters::ppm4;

  for(std::size_t axis = 0; axis < D; axis++) {
    // clang-format off
    flecsi::execute<tasks::hydro::reconstruct<D, limiter>, flecsi::default_accelerator>(axis,
      s.m, s.mass_density(s.m), s.velocity(s.m), s.pressure(s.m),
      s.radiation_energy_density(s.m),
      s.rTail(s.m), s.rHead(s.m), s.uTail(s.m), s.uHead(s.m),
      s.pTail(s.m), s.pHead(s.m), s.EradTail(s.m), s.EradHead(s.m),
      s.ruTail(s.m), s.ruHead(s.m), s.rETail(s.m), s.rEHead(s.m), gamma(s.gt));

    if constexpr(Stage == time_stepper::rk_stage::First) {
      flecsi::execute<tasks::hydro::compute_interface_fluxes<D>, flecsi::default_accelerator>(axis, s.m,
        s.rTail(s.m), s.rHead(s.m), s.uTail(s.m), s.uHead(s.m),
        s.pTail(s.m), s.pHead(s.m), s.EradTail(s.m), s.EradHead(s.m),
        s.ruTail(s.m), s.ruHead(s.m),
        s.rETail(s.m), s.rEHead(s.m),
        s.rF(s.m), s.ruF(s.m), s.rEF(s.m), s.EradF(s.m),
        s.dt_mass_density(s.m),
        s.dt_momentum_density(s.m),
        s.dt_total_energy_density(s.m),
        s.dt_radiation_energy_density(s.m),
        gamma(s.gt));
    }
    else if constexpr(Stage == time_stepper::rk_stage::Second) {
      flecsi::execute<tasks::hydro::compute_interface_fluxes<D>, flecsi::default_accelerator>(axis, s.m,
        s.rTail(s.m), s.rHead(s.m), s.uTail(s.m),s.uHead(s.m),
        s.pTail(s.m), s.pHead(s.m), s.EradTail(s.m), s.EradHead(s.m),
        s.ruTail(s.m), s.ruHead(s.m), s.rETail(s.m), s.rEHead(s.m),
        s.rF(s.m), s.ruF(s.m), s.rEF(s.m), s.EradF(s.m),
        s.dt_mass_density_2(s.m),
        s.dt_momentum_density_2(s.m),
        s.dt_total_energy_density_2(s.m),
        s.dt_radiation_energy_density_2(s.m),
        gamma(s.gt));
    }
    // clang-format on
  }
} // fluxes_terms

// --------------------------------------------------------------------
//              IMPLICIT PARTS OF EVOLUTION
// --------------------------------------------------------------------

// radiative heating-cooling part
template<std::size_t D, time_stepper::rk_stage Stage>
void
implicit_source_terms(control_policy<state, D> & cp) {

#ifndef DISABLE_RADIATION

  auto & s = cp.state();

  // Compute and store (de_dt, dE_dt)_implicit from matter-radiation coupling.
  if constexpr(Stage == time_stepper::rk_stage::First) {
    flecsi::execute<task::rad_root::update_energy_density<D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.velocity(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      kappa(s.gt),
      particle_mass(s.gt),
      s.dt_weighted(s.gt),
      s.dt_total_energy_density_implicit(s.m),
      s.dt_radiation_energy_density_implicit(s.m));
  }
  else if constexpr(Stage == time_stepper::rk_stage::Second) {
    flecsi::execute<task::rad_root::update_energy_density<D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.velocity(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      kappa(s.gt),
      particle_mass(s.gt),
      s.dt_weighted(s.gt),
      s.dt_total_energy_density_implicit_2(s.m),
      s.dt_radiation_energy_density_implicit_2(s.m));
  }

#endif
}

// --------------------------------------------------------------------
//  Update conserved variables with computed (dU_dt)_explicit and
//  (dU_dt)_implicit. Then perform primitive recovery and apply boundary
//  conditions
// --------------------------------------------------------------------
template<std::size_t D, time_stepper::rk_stage Stage>
void
update_variables(control_policy<state, D> & cp) {
  auto & s = cp.state();

  flecsi::execute<tasks::update_u<D, Stage>, flecsi::default_accelerator>(
    s.dt(s.gt),
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
    s.dt_radiation_energy_density(s.m),
    s.dt_total_energy_density_implicit(s.m),
    s.dt_radiation_energy_density_implicit(s.m),
    //
    s.dt_mass_density_2(s.m),
    s.dt_momentum_density_2(s.m),
    s.dt_total_energy_density_2(s.m),
    s.dt_radiation_energy_density_2(s.m),
    s.dt_total_energy_density_implicit_2(s.m),
    s.dt_radiation_energy_density_implicit_2(s.m),
    //
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.radiation_energy_density(s.m));

  // Perform primitive recovery
  flecsi::execute<tasks::hydro::conservative_to_primitive<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    gamma(s.gt));

  // Update boundary cells
  flecsi::execute<tasks::apply_boundaries<D>, flecsi::default_accelerator>(s.m,
    s.bmap(s.gt),
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.radiation_energy_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m));

} // radiation_advance

// --------------------------------------------------------------------
//  Compute max characteristic speeds, and determine dt_min() for the next time
//  step.
// --------------------------------------------------------------------
template<std::size_t D>
void
update_time_step_size(control_policy<state, D> & cp) {
  auto & s = cp.state();

  auto lmax_f =
    flecsi::execute<tasks::hydro::update_max_characteristic_speed<D>,
      flecsi::default_accelerator>(
      s.m, s.mass_density(s.m), s.velocity(s.m), s.pressure(s.m), gamma(s.gt));

  s.dtmin_ =
    flecsi::reduce<tasks::hydro::update_dtmin<D>, flecsi::exec::fold::min>(
      s.m, lmax_f);

#ifdef HARD_ENABLE_LEGION_TRACING
  cp.guard.reset();
#endif
}

} // namespace hard::actions
