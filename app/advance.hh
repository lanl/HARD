#ifndef HARD_ADVANCE_HH
#define HARD_ADVANCE_HH

#include "actions.hh"
#include "numerical_algorithms/time_stepper.hh"
#include "rad.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/rad.hh"
#include "tasks/time_derivative.hh"

#include <flecsi/execution.hh>

namespace hard::actions {

template<std::size_t D, time_stepper::rk_stage Stage>
void
radiation_advance(control_policy<state, D> & cp) {
  using namespace flecsi;
  auto & s = cp.state();

#ifndef DISABLE_RADIATION

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

  // After solving the elliptic part, store dU_dt for implicit part and
  // recover the intermediate state
  if constexpr(Stage == time_stepper::rk_stage::First) {
    flecsi::execute<task::rad::store_du_dt_implicit_from_diffusion<D>,
      flecsi::default_accelerator>(s.m,
      s.Esf(s.m),
      s.Ef(s.m),
      s.dt_weighted(s.gt),
      s.dt_radiation_energy_density_implicit(s.m));

    flecsi::execute<tasks::compute_u_after_implicit_solve<D>,
      flecsi::default_accelerator>(s.m,
      s.dt_weighted(s.gt),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.dt_total_energy_density_implicit(s.m),
      s.dt_radiation_energy_density_implicit(s.m));
  }
  else if constexpr(Stage == time_stepper::rk_stage::Second) {
    flecsi::execute<task::rad::store_du_dt_implicit_from_diffusion<D>,
      flecsi::default_accelerator>(s.m,
      s.Esf(s.m),
      s.Ef(s.m),
      s.dt_weighted(s.gt),
      s.dt_radiation_energy_density_implicit_2(s.m));

    flecsi::execute<tasks::compute_u_after_implicit_solve<D>,
      flecsi::default_accelerator>(s.m,
      s.dt_weighted(s.gt),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.dt_total_energy_density_implicit_2(s.m),
      s.dt_radiation_energy_density_implicit_2(s.m));
  }

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
#endif

  // and also update boundary cells
  flecsi::execute<tasks::apply_boundaries<D>, flecsi::default_accelerator>(s.m,
    s.bmap(s.gt),
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.radiation_energy_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m));

} // radiation_advance

// ---------------------------------------------------------------------------------------
//              Action lists and dependencies
// ---------------------------------------------------------------------------------------

using namespace hard::time_stepper;

//
//  < Second order IMEX-RK time stepping >
//

//  Reference: IMEX-SSP2(2,2,2) scheme, Table 2 of Pareschi & Russo 2005
//  (https://arxiv.org/abs/1009.2757)
//
//  0. Start with the current state of evolved variables, U^n
//  1. Implicit solve for U^1 = U^n + g I(U^1).
//     Store E(U^1) and I(U^1)
//  2. Implicit solve for U^2 = U^n + E(U^1) + (1-2g) * I(U^1) + g * I(U^2)
//     Store E(U^2) and I(U^2).
//  3. U^{n+1} = U^n + 0.5 * (E(U^1) + E(U^2) + I(U^1) + I(U^2))
//

inline control<state, 1>::action<initialize_time_derivative<1>, cp::advance>
  initialize_time_derivative_1d;
inline control<state, 2>::action<initialize_time_derivative<2>, cp::advance>
  initialize_time_derivative_2d;
inline control<state, 3>::action<initialize_time_derivative<3>, cp::advance>
  initialize_time_derivative_3d;

// --------------------------------------------------------------------
//              RK substep 1
// --------------------------------------------------------------------

// Implicit solve on U^n to get U^1
//  a) Heating
inline control<state, 1>::action<implicit_source_terms<1, rk_stage::First>,
  cp::advance>
  implicit_source_terms_first_1d;
inline control<state, 2>::action<implicit_source_terms<2, rk_stage::First>,
  cp::advance>
  implicit_source_terms_first_2d;
inline control<state, 3>::action<implicit_source_terms<3, rk_stage::First>,
  cp::advance>
  implicit_source_terms_first_3d;
inline const auto dep_implicit_source_terms_first_1d =
  implicit_source_terms_first_1d.add(initialize_time_derivative_1d);
inline const auto dep_implicit_source_terms_first_2d =
  implicit_source_terms_first_2d.add(initialize_time_derivative_2d);
inline const auto dep_implicit_source_terms_first_3d =
  implicit_source_terms_first_3d.add(initialize_time_derivative_3d);
//  b) Diffusion & recover U^1. This also stores I(U^1)
inline control<state, 1>::action<radiation_advance<1, rk_stage::First>,
  cp::advance>
  radiation_advance_first_1d;
inline control<state, 2>::action<radiation_advance<2, rk_stage::First>,
  cp::advance>
  radiation_advance_first_2d;
inline control<state, 3>::action<radiation_advance<3, rk_stage::First>,
  cp::advance>
  radiation_advance_first_3d;
inline const auto dep_radiation_advance_first_1d =
  radiation_advance_first_1d.add(implicit_source_terms_first_1d);
inline const auto dep_radiation_advance_first_2d =
  radiation_advance_first_2d.add(implicit_source_terms_first_2d);
inline const auto dep_radiation_advance_first_3d =
  radiation_advance_first_3d.add(implicit_source_terms_first_3d);

// Compute and store explicit part E(U^1)
//  a) flux part
inline control<state, 1>::action<explicit_source_terms<1, rk_stage::First>,
  cp::advance>
  explicit_source_terms_first_1d;
inline control<state, 2>::action<explicit_source_terms<2, rk_stage::First>,
  cp::advance>
  explicit_source_terms_first_2d;
inline control<state, 3>::action<explicit_source_terms<3, rk_stage::First>,
  cp::advance>
  explicit_source_terms_first_3d;
inline const auto dep_explicit_source_terms_first_1d =
  explicit_source_terms_first_1d.add(radiation_advance_first_1d);
inline const auto dep_explicit_source_terms_first_2d =
  explicit_source_terms_first_2d.add(radiation_advance_first_2d);
inline const auto dep_explicit_source_terms_first_3d =
  explicit_source_terms_first_3d.add(radiation_advance_first_3d);
//  b) explicit source part
inline control<state, 1>::action<fluxes_terms<1, rk_stage::First>, cp::advance>
  fluxes_terms_first_1d;
inline control<state, 2>::action<fluxes_terms<2, rk_stage::First>, cp::advance>
  fluxes_terms_first_2d;
inline control<state, 3>::action<fluxes_terms<3, rk_stage::First>, cp::advance>
  fluxes_terms_first_3d;
inline const auto dep_fluxes_terms_first_1d =
  fluxes_terms_first_1d.add(explicit_source_terms_first_1d);
inline const auto dep_fluxes_terms_first_2d =
  fluxes_terms_first_2d.add(explicit_source_terms_first_2d);
inline const auto dep_fluxes_terms_first_3d =
  fluxes_terms_first_3d.add(explicit_source_terms_first_3d);

// --------------------------------------------------------------------
//              RK substep 2
// --------------------------------------------------------------------

// Compute U^* = U^n + E(U^1) + (1-2g) * I(U^1)
inline control<state, 1>::action<update_variables<1, rk_stage::Second>,
  cp::advance>
  update_vars_second_1d;
inline control<state, 2>::action<update_variables<2, rk_stage::Second>,
  cp::advance>
  update_vars_second_2d;
inline control<state, 3>::action<update_variables<3, rk_stage::Second>,
  cp::advance>
  update_vars_second_3d;
inline const auto dep_update_vars_second_1d =
  update_vars_second_1d.add(fluxes_terms_first_1d);
inline const auto dep_update_vars_second_2d =
  update_vars_second_2d.add(fluxes_terms_first_2d);
inline const auto dep_update_vars_second_3d =
  update_vars_second_3d.add(fluxes_terms_first_3d);

// Implicit solve on U^* to get U^2
//  a) Heating
inline control<state, 1>::action<implicit_source_terms<1, rk_stage::Second>,
  cp::advance>
  implicit_source_terms_second_1d;
inline control<state, 2>::action<implicit_source_terms<2, rk_stage::Second>,
  cp::advance>
  implicit_source_terms_second_2d;
inline control<state, 3>::action<implicit_source_terms<3, rk_stage::Second>,
  cp::advance>
  implicit_source_terms_second_3d;
inline const auto dep_implicit_source_terms_second_1d =
  implicit_source_terms_second_1d.add(update_vars_second_1d);
inline const auto dep_implicit_source_terms_second_2d =
  implicit_source_terms_second_2d.add(update_vars_second_2d);
inline const auto dep_implicit_source_terms_second_3d =
  implicit_source_terms_second_3d.add(update_vars_second_3d);
//  b) Diffusion & recover U^2. This also stores I(U^2)
inline control<state, 1>::action<radiation_advance<1, rk_stage::Second>,
  cp::advance>
  radiation_advance_second_1d;
inline control<state, 2>::action<radiation_advance<2, rk_stage::Second>,
  cp::advance>
  radiation_advance_second_2d;
inline control<state, 3>::action<radiation_advance<3, rk_stage::Second>,
  cp::advance>
  radiation_advance_second_3d;
inline const auto dep_radiation_advance_second_1d =
  radiation_advance_second_1d.add(implicit_source_terms_second_1d);
inline const auto dep_radiation_advance_second_2d =
  radiation_advance_second_2d.add(implicit_source_terms_second_2d);
inline const auto dep_radiation_advance_second_3d =
  radiation_advance_second_3d.add(implicit_source_terms_second_3d);

// Compute and store explicit part E(U^2)
//  a) flux part
inline control<state, 1>::action<explicit_source_terms<1, rk_stage::Second>,
  cp::advance>
  explicit_source_terms_second_1d;
inline control<state, 2>::action<explicit_source_terms<2, rk_stage::Second>,
  cp::advance>
  explicit_source_terms_second_2d;
inline control<state, 3>::action<explicit_source_terms<3, rk_stage::Second>,
  cp::advance>
  explicit_source_terms_second_3d;
inline const auto dep_explicit_source_terms_second_1d =
  explicit_source_terms_second_1d.add(radiation_advance_second_1d);
inline const auto dep_explicit_source_terms_second_2d =
  explicit_source_terms_second_2d.add(radiation_advance_second_2d);
inline const auto dep_explicit_source_terms_second_3d =
  explicit_source_terms_second_3d.add(radiation_advance_second_3d);
//  b) explicit source part
inline control<state, 1>::action<fluxes_terms<1, rk_stage::Second>, cp::advance>
  fluxes_terms_second_1d;
inline control<state, 2>::action<fluxes_terms<2, rk_stage::Second>, cp::advance>
  fluxes_terms_second_2d;
inline control<state, 3>::action<fluxes_terms<3, rk_stage::Second>, cp::advance>
  fluxes_terms_second_3d;
inline const auto dep_fluxes_terms_second_1d =
  fluxes_terms_second_1d.add(explicit_source_terms_second_1d);
inline const auto dep_fluxes_terms_second_2d =
  fluxes_terms_second_2d.add(explicit_source_terms_second_2d);
inline const auto dep_fluxes_terms_second_3d =
  fluxes_terms_second_3d.add(explicit_source_terms_second_3d);

// --------------------------------------------------------------------
//              RK substep 3
// --------------------------------------------------------------------

// Perform the full time step U^n -> U^n+1
inline control<state, 1>::action<update_variables<1, rk_stage::Update>,
  cp::advance>
  perform_time_step_1d;
inline control<state, 2>::action<update_variables<2, rk_stage::Update>,
  cp::advance>
  perform_time_step_2d;
inline control<state, 3>::action<update_variables<3, rk_stage::Update>,
  cp::advance>
  perform_time_step_3d;
inline const auto dep_perform_time_step_1d =
  perform_time_step_1d.add(fluxes_terms_second_1d);
inline const auto dep_perform_time_step_2d =
  perform_time_step_2d.add(fluxes_terms_second_2d);
inline const auto dep_perform_time_step_3d =
  perform_time_step_3d.add(fluxes_terms_second_3d);

// Update the CFL limit
inline control<state, 1>::action<update_time_step_size<1>, cp::advance>
  update_time_step_size_1d;
inline control<state, 2>::action<update_time_step_size<2>, cp::advance>
  update_time_step_size_2d;
inline control<state, 3>::action<update_time_step_size<3>, cp::advance>
  update_time_step_size_3d;
inline const auto dep_update_time_step_size_1d =
  update_time_step_size_1d.add(perform_time_step_1d);
inline const auto dep_update_time_step_size_2d =
  update_time_step_size_2d.add(perform_time_step_2d);
inline const auto dep_update_time_step_size_3d =
  update_time_step_size_3d.add(perform_time_step_3d);

} // namespace hard::actions

#endif // HARD_ADVANCE_HH
