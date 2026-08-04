#ifndef HARD_APPS_ACTIONS_RK_STAGES_HH
#define HARD_APPS_ACTIONS_RK_STAGES_HH

#include "state.hh"

#include "gfm/tasks/cons2prim.hh"
#include "gfm/tasks/extrapolation.hh"
#include "gfm/tasks/interface_fluxes_gfm.hh"
#include "gfm/tasks/reconstruct.hh"
#include "gfm/tasks/utils.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

#include "gfm/tasks/compute_sound_speed.hh"

#include <algorithm>

namespace hard {

template<std::size_t D>
using field_double =
  flecsi::field<double>::Reference<spec::mesh<D>, spec::is::cells>;
template<std::size_t D>
using face_pair = typename hard::faces<D>::ref_pair;

//----------------------------------------------------------------------------//
// Shared building blocks for the multi-material hydro RK stages.
//----------------------------------------------------------------------------//

// Perform primitive recovery from the conservative variables.
template<std::size_t D>
inline void
perform_cons2prim(state<D> & s, flecsi::scheduler & sc) {
  annotation_rguard<tasks::gfm::conservative_to_primitive_mm_task>(),
    sc.execute<tasks::gfm::conservative_to_primitive_mm<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.prim.velocity(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.prim.sound_speed(*s.m),
      s.ls.material(*s.m),
      s.eos_mm);
}

// Update the INTEGER material field and the DOUBLE psi field in ghost cells.
template<std::size_t D>
inline void
apply_material_boundaries(state<D> & s, flecsi::scheduler & sc) {
  // Update INTEGER material field in ghost cells before hydro/primitive
  // variables.
  annotation_rguard<tasks::apply_boundaries_scalar_task>(),
    sc.execute<tasks::apply_boundaries_scalar<D, true, int>>(flecsi::exec::on,
      *s.m,
      s.icst.bmap(*s.gt),
      std::vector{s.ls.material(*s.m)}); // Only int fields

  // Update DOUBLE psi field in ghost cells
  annotation_rguard<tasks::apply_boundaries_scalar_task>(),
    sc.execute<tasks::apply_boundaries_scalar<D, true, double>>(
      flecsi::exec::on,
      *s.m,
      s.icst.bmap(*s.gt),
      std::vector{s.ls.psi(*s.m)}); // Only double scalar fields
}

// Update material/psi ghost cells followed by the primitive variables in ghost
// cells for the supplied set of scalar fields.
template<std::size_t D>
inline void
apply_scalar_boundaries(state<D> & s,
  flecsi::scheduler & sc,
  std::vector<field_double<D>> scalar_v) {
  apply_material_boundaries(s, sc);

  // Update boundary cells
  annotation_rguard<tasks::apply_boundaries_task>(),
    sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
      *s.m,
      s.icst.bmap(*s.gt),
      scalar_v,
      std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});
}

// Recompute sound speed from EOS in ghost cells.
template<std::size_t D>
inline void
recompute_sound_speed(state<D> & s, flecsi::scheduler & sc) {
  annotation_rguard<hard::tasks::gfm::compute_sound_speed_from_eos_task>(),
    sc.execute<hard::tasks::gfm::compute_sound_speed_from_eos<D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.ls.material(*s.m),
      s.prim.sound_speed(*s.m),
      s.eos_mm);
}

// Aslam JCP 2003 constant extrapolation routine used to fill the ghost stage
// for a single source material (0 or 1 in the current two-material setup).
// The result is written into q_ext_out.
template<std::size_t D>
inline void
extrapolate_all_materials(state<D> & s, flecsi::scheduler & sc) {

  // 1) Extrapolate every material into its own key.
  for(std::size_t mat = 0; mat < s.n_materials(); ++mat) {

    auto q_ext_out =
      (mat == 0) ? s.extp.q_ext_plus(*s.m) : s.extp.q_ext_minus(*s.m);

    sc.execute<tasks::extrapolation::prepare_extp_for_material<D>>(
      flecsi::exec::on,
      *s.m,
      s.ls.psi(*s.m),
      s.ls.material(*s.m),
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.extp.extp_psi(*s.m),
      q_ext_out,
      s.extp.q_scalar_known(*s.m),
      static_cast<int>(mat));

    // Compute unit normal of level set gradient
    // dp = grad(psi)
    annotation_rguard<tasks::extrapolation::gradient_centered_task>(),
      sc.execute<tasks::extrapolation::gradient_centered<D>>(
        flecsi::exec::on, *s.m, s.extp.extp_psi(*s.m), s.extp.dp(*s.m));

    // normal = dp/|dp| and normal faces
    annotation_rguard<tasks::extrapolation::unit_cell_normals_task>(),
      sc.execute<tasks::extrapolation::unit_cell_normals<D>>(
        flecsi::exec::on, *s.m, s.extp.dp(*s.m), s.extp.normal(*s.m));

    for(int n = 0; n < s.extp_steps; ++n) {
      // Stage 1: rhs(q^n) -> rhs
      annotation_rguard<
        tasks::extrapolation::rhs_constant_extrapolation_task>(),
        sc.execute<tasks::extrapolation::rhs_constant_extrapolation<D>>(
          flecsi::exec::on,
          *s.m,
          q_ext_out,
          s.extp.extp_psi(*s.m),
          s.extp.normal(*s.m),
          s.extp.q_rhs(*s.m));

      // q1 = q^n + dt * rhs(q^n)
      annotation_rguard<tasks::extrapolation::advance_q_rk2_s1_cfl_task>(),
        sc.execute<tasks::extrapolation::advance_q_rk2_s1_cfl<D>>(
          flecsi::exec::on,
          *s.m,
          s.extp.extp_psi(*s.m),
          s.extp.normal(*s.m),
          s.extp.q_scalar_known(*s.m),
          q_ext_out,
          s.extp.q_rhs(*s.m),
          s.extp.q_scalar1(*s.m),
          tasks::extrapolation::extp_cfl);

      // Stage 2: rhs(q1) -> rhs1
      annotation_rguard<
        tasks::extrapolation::rhs_constant_extrapolation_task>(),
        sc.execute<tasks::extrapolation::rhs_constant_extrapolation<D>>(
          flecsi::exec::on,
          *s.m,
          s.extp.q_scalar1(*s.m),
          s.extp.extp_psi(*s.m),
          s.extp.normal(*s.m),
          s.extp.q_rhs1(*s.m));

      // q^{n+1} = 0.5*( q^n + q1 + dt*rhs(q1) )
      annotation_rguard<tasks::extrapolation::advance_q_rk2_s2_cfl_task>(),
        sc.execute<tasks::extrapolation::advance_q_rk2_s2_cfl<D>>(
          flecsi::exec::on,
          *s.m,
          s.extp.extp_psi(*s.m),
          s.extp.normal(*s.m),
          s.extp.q_scalar_known(*s.m),
          q_ext_out,
          s.extp.q_rhs1(*s.m),
          s.extp.q_scalar1(*s.m),
          tasks::extrapolation::extp_cfl); // dtau in original work

      // Copy extrapolated q_ext_out into q_ext_tmp to check convergence
      if(n == 0 || n == (s.extp_steps - 2)) {
        annotation_rguard<tasks::extrapolation::copy_extp_q_task>(),
          sc.execute<tasks::extrapolation::copy_extp_q<D>>(
            flecsi::exec::on, q_ext_out, s.extp.q_scalar_tmp(*s.m));
      }
    }
    // Post-check the pseudo-time extrapolation using the active SciPy-style
    // scaled update metric in the first extrapolated layer next to the
    // interface.
    auto max_scaled_update_band_f =
      (annotation_rguard<
         tasks::extrapolation::extp_interface_convergence_update_task>(),
        sc.execute<tasks::extrapolation::extp_interface_convergence_update<D>>(
          flecsi::exec::on,
          *s.m,
          s.extp.q_scalar_tmp(*s.m),
          q_ext_out,
          s.extp.extp_psi(*s.m)));

    // check here if the extrapolation has converged
    annotation_rguard<tasks::extrapolation::extp_interface_convergence_task>(),
      sc.execute<tasks::extrapolation::extp_interface_convergence<D>>(
        flecsi::exec::on,
        *s.m,
        max_scaled_update_band_f,
        mat,
        s.extp_steps,
        tasks::extrapolation::extp_residual_tol);
  }
}

// Per-axis loop shared by both RK stages: extrapolation (ghost fill),
// reconstruction and interface flux computation. Note that the interface
// fluxes always accumulate into rk_dt1 (K1)
template<std::size_t D, class F>
inline void
advance_rk_stage_axis_loop(state<D> & s,
  flecsi::scheduler & sc,
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> & v_t,
  F f) {

  using limiter = spec::limiters::weno5z;

  for(std::size_t axis = 0; axis < D; axis++) {
    {
      auto scalar_v = std::vector{
        std::make_tuple(s.cons.hydro.mass_density(*s.m), s.f.hydro.r_face(s.m)),
        std::make_tuple(
          s.prim.specific_internal_energy(*s.m), s.f.hydro.e_face(s.m)),
        std::make_tuple(s.prim.sound_speed(*s.m), s.f.hydro.c_face(s.m)),
        std::make_tuple(s.prim.pressure(*s.m), s.f.hydro.p_face(s.m))};

      for(auto & v : v_t) {
        scalar_v.push_back(std::make_tuple(std::get<0>(v), std::get<2>(v)));
      }

      // clang-format off
      annotation_rguard<tasks::hydro::reconstruct_primitives_task>(),
        sc.execute<tasks::hydro::reconstruct_primitives<D, limiter>>(
          flecsi::exec::on,
          axis,
          *s.m,
          s.ls.material(*s.m),
          scalar_v,
          std::vector{
            std::make_tuple(
              s.prim.velocity(*s.m),
              s.f.hydro.u_face(s.m))});
      // clang-format on
    }

    annotation_rguard<tasks::hydro::reconstruct_conserved_task>(),
      sc.execute<tasks::hydro::reconstruct_conserved<D>>(flecsi::exec::on,
        *s.m,
        s.f.hydro.r_face(s.m),
        s.f.hydro.u_face(s.m),
        s.f.hydro.e_face(s.m),
        s.f.hydro.ru_face(s.m),
        s.f.hydro.re_face(s.m));

    // Calculate K1 and save it to dt_U
    // Runtime dispatch based on GFM method
    if(s.extp.gfm_method == state<D>::extrapolation::original_gfm) {
      annotation_rguard<tasks::gfm::compute_interface_fluxes_gfm_task>(),
        sc.execute<tasks::gfm::compute_interface_fluxes_gfm<D,
          state<D>::extrapolation::original_gfm>>(flecsi::exec::on,
          axis,
          *s.m,
          s.f.hydro.r_face(s.m),
          s.f.hydro.u_face(s.m),
          s.f.hydro.p_face(s.m),
          s.f.hydro.c_face(s.m),
          s.f.hydro.ru_face(s.m),
          s.f.hydro.re_face(s.m),
          // GFM
          s.ls.material(*s.m),
          s.extp.normal(*s.m),
          s.extp.q_ext_plus(*s.m),
          s.extp.q_ext_minus(*s.m),
          // Riemann Fluxes
          s.rf.hydro.r_f(*s.m),
          s.rf.hydro.ru_f(*s.m),
          s.rf.hydro.re_f(*s.m),
          s.rk_dt1(s.m),
          s.icst.gravity_acc(*s.gt),
          s.eos_mm);
    }
    else {
      annotation_rguard<tasks::gfm::compute_interface_fluxes_gfm_task>(),
        sc.execute<tasks::gfm::compute_interface_fluxes_gfm<D,
          state<D>::extrapolation::riemann_gfm>>(flecsi::exec::on,
          axis,
          *s.m,
          s.f.hydro.r_face(s.m),
          s.f.hydro.u_face(s.m),
          s.f.hydro.p_face(s.m),
          s.f.hydro.c_face(s.m),
          s.f.hydro.ru_face(s.m),
          s.f.hydro.re_face(s.m),
          // GFM
          s.ls.material(*s.m),
          s.extp.normal(*s.m),
          s.extp.q_ext_plus(*s.m),
          s.extp.q_ext_minus(*s.m),
          // Riemann Fluxes
          s.rf.hydro.r_f(*s.m),
          s.rf.hydro.ru_f(*s.m),
          s.rf.hydro.re_f(*s.m),
          s.rk_dt1(s.m),
          s.icst.gravity_acc(*s.gt),
          s.eos_mm);
    }

    // Apply function if present
    if constexpr(!std::is_same_v<F, std::nullptr_t>)
      f(axis);
  }
}

// Common driver for an explicit RK stage. rk_dt selects the time-derivative
// storage used for the explicit source term and the initial conservative
// update (rk_dt1 for stage 1, rk_dt2 for stage 2). do_sound_speed toggles the
// EOS sound-speed recomputation in ghost cells (enabled for stage 1 only).
template<std::size_t D, class F = std::nullptr_t>
inline void
advance_rk_stage_common(state<D> & s,
  flecsi::scheduler & sc,
  RK<D> & rk_dt,
  bool do_sound_speed,
  // conservative, dt, face
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> v_t =
    {},
  F f = nullptr) {

  // Explicit source term (gravity) update for RT case in hydro file
  annotation_rguard<tasks::external_source_task>(),
    sc.execute<tasks::external_source<D>>(flecsi::exec::on,
      *s.m,
      s.prim.velocity(*s.m),
      s.src_t.hydro.gravity_force(*s.m),
      // time-derivatives
      rk_dt.momentum_energy_density()(*s.m),
      rk_dt.total_energy_density()(*s.m));

  {
    auto scalar_v = std::vector{std::make_tuple(rk_dt.mass_density()(*s.m),
                                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(rk_dt.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(std::make_tuple(std::get<1>(v), std::get<0>(v)));
    }

    // We need update_u here before we compute fluxes in the presence of
    // hydro::explict_source_update with body forces. See Moens'21 Eq. 24-26
    // clang-format off
    annotation_rguard<tasks::hydro::update_u_task>(),
      sc.execute<tasks::hydro::update_u<D>>(flecsi::exec::on,
        s.dt(*s.gt),
        scalar_v,
        std::vector{
          std::make_tuple(
            rk_dt.momentum_energy_density()(*s.m),
            s.cons.hydro.momentum_density(*s.m))});
    // clang-format on
  }

  // Perform primitive recovery
  perform_cons2prim(s, sc);

  // Update material/psi and primitive variables in ghost cells BEFORE
  // reconstruction.
  apply_scalar_boundaries(s,
    sc,
    std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      // s.prim.sound_speed(*s.m),
      s.cons.hydro.total_energy_density(*s.m)});

  // Recompute sound speed from EOS in ghost cells (stage 1 only).
  if(do_sound_speed)
    recompute_sound_speed(s, sc);

  /* Start extrapolation routine to fill ghost stage
   * This is based on Aslam JCP 2003, constant extrapolation */

  // Extrapolation for ghost fill
  extrapolate_all_materials(s, sc);
  /* End of extrapolation*/

  advance_rk_stage_axis_loop(s, sc, v_t, f);
}

} // namespace hard

#endif // HARD_APPS_ACTIONS_RK_STAGES_HH
