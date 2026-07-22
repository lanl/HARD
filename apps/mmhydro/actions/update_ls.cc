#include "state.hh"

#include "gfm/tasks/ls_core.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

#ifdef IN_TEST
#include "tests_level_set_utils.hh"
#endif

namespace hard {

struct update_ls_ssprk3_region : annotation_action_region {
  inline static const std::string name{"update_ls_ssprk3"};
};

using namespace hard::tasks::level_set;

template<std::size_t D>
void
ssprk3_step(
#ifndef IN_TEST
  control_policy<hard::state, D> & cp
#else
  flecsi::scheduler & sc,
  tests::level_set_utils::state<D> & s
#endif
) {

  annotation_rguard<update_ls_ssprk3_region> guard;

#ifndef IN_TEST
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();
#endif

#ifdef IN_TEST
  sc.execute<tests::level_set_utils::copy_test_vel_stage<D, 1>>(
    flecsi::exec::on,
    *s.m,
    s.velocity(*s.m),
    s.t(*s.gt),
    s.dt(*s.gt),
    s.period(*s.gt));
#endif // IN_TEST

  // 0) Freeze material for this entire time step: Lfrozen = L^n
  annotation_rguard<freeze_material_task>(),
    sc.execute<freeze_material<D>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.material(*s.m),
      s.frozen_material(*s.m));
#else
      s.ls.material(*s.m),
      s.ls.frozen_material(*s.m));
#endif

  // Stage 1: k = L(psi^n, L^n); psi1 = psi^n + dt*k
  annotation_rguard<rhs_upwind_task>(),
    sc.execute<rhs_upwind<D>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi(*s.m),
      s.velocity(*s.m),
      s.frozen_material(*s.m),
      s.rhs(*s.m));
#else
      s.ls.psi(*s.m),
      s.prim.velocity(*s.m),
      s.ls.frozen_material(*s.m),
      s.ls.rhs(*s.m));
#endif

  annotation_rguard<advance_psi_s1_task>(),
    sc.execute<advance_psi_s1<D, 1>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi(*s.m), // psi^n
      s.rhs(*s.m), // k(psi^n)
      s.psi1(*s.m), // psi1
#else
      s.ls.psi(*s.m), // psi^n
      s.ls.rhs(*s.m), // k(psi^n)
      s.ls.psi1(*s.m), // psi1
#endif
      s.dt(*s.gt));

#ifdef IN_TEST
  sc.execute<tests::level_set_utils::copy_test_vel_stage<D, 2>>(
    flecsi::exec::on,
    *s.m,
    s.velocity(*s.m),
    s.t(*s.gt),
    s.dt(*s.gt),
    s.period(*s.gt));
#endif // IN_TEST

  // Stage 2: k = L(psi1, L^n); psi2 = 0.75*psi^n + 0.25*(psi1 + dt*k)
  annotation_rguard<rhs_upwind_task>(),
    sc.execute<rhs_upwind<D>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi1(*s.m),
      s.velocity(*s.m),
      s.frozen_material(*s.m),
      s.rhs(*s.m));
#else
      s.ls.psi1(*s.m),
      s.prim.velocity(*s.m),
      s.ls.frozen_material(*s.m),
      s.ls.rhs(*s.m));
#endif

  annotation_rguard<advance_psi_task>(),
    sc.execute<advance_psi<D, 2>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi(*s.m), // psi^n
      s.rhs(*s.m), // k(psi1)
      s.psi1(*s.m), // psi_prev = psi1
      s.psi2(*s.m), // psi2
#else
      s.ls.psi(*s.m), // psi^n
      s.ls.rhs(*s.m), // k(psi1)
      s.ls.psi1(*s.m), // psi_prev = psi1
      s.ls.psi2(*s.m), // psi2
#endif
      s.dt(*s.gt));

#ifdef IN_TEST
  sc.execute<tests::level_set_utils::copy_test_vel_stage<D, 3>>(
    flecsi::exec::on,
    *s.m,
    s.velocity(*s.m),
    s.t(*s.gt),
    s.dt(*s.gt),
    s.period(*s.gt));
#endif // IN_TEST

  // Stage 3: k = L(psi2, L^n); psi3 = 1/3*psi^n + 2/3*(psi2 + dt*k)
  annotation_rguard<rhs_upwind_task>(),
    sc.execute<rhs_upwind<D>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi2(*s.m),
      s.velocity(*s.m),
      s.frozen_material(*s.m),
      s.rhs(*s.m));
#else
      s.ls.psi2(*s.m),
      s.prim.velocity(*s.m),
      s.ls.frozen_material(*s.m),
      s.ls.rhs(*s.m));
#endif

  annotation_rguard<advance_psi_task>(),
    sc.execute<advance_psi<D, 3>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi(*s.m), // psi^n
      s.rhs(*s.m), // k(psi2)
      s.psi2(*s.m), // psi_prev = psi2
      s.psi3(*s.m), // psi3
#else
      s.ls.psi(*s.m), // psi^n
      s.ls.rhs(*s.m), // k(psi2)
      s.ls.psi2(*s.m), // psi_prev = psi2
      s.ls.psi3(*s.m), // psi3
#endif
      s.dt(*s.gt));

  // Algorithm 1 post-process: (psi3, L^n) -> (psi^{n+1}, L^{n+1})
  annotation_rguard<determine_psi_task>(),
    sc.execute<determine_psi<D>>(flecsi::exec::on,
      *s.m,
#ifdef IN_TEST
      s.psi3(*s.m), // frozen_psi = psi3
      s.psi(*s.m), // output psi^{n+1}
      s.frozen_material(*s.m), // L^n
      s.material(*s.m)); // output L^{n+1}
#else
      s.ls.psi3(*s.m), // frozen_psi = psi3
      s.ls.psi(*s.m), // output psi^{n+1}
      s.ls.frozen_material(*s.m), // L^n
      s.ls.material(*s.m)); // output L^{n+1}
#endif

// ------ ADDED NEW 1 ------- START
#ifndef IN_TEST
  {
    annotation_rguard<promote_changed_material_cells_task>(),
      sc.execute<promote_changed_material_cells<D>>(flecsi::exec::on,
        *s.m,
        s.ls.frozen_material(*s.m), // old material L^n
        s.ls.material(*s.m), // new material L^{n+1}
        s.extp.q_ext_plus(*s.m), // material 0 extrapolated state
        s.extp.q_ext_minus(*s.m), // material 1 extrapolated state
        s.cons.hydro.mass_density(*s.m),
        s.cons.hydro.momentum_density(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.prim.velocity(*s.m),
        s.prim.pressure(*s.m),
        s.prim.specific_internal_energy(*s.m),
        s.prim.sound_speed(*s.m),
        s.eos_mm);

    {
      auto scalar_v = std::vector{s.cons.hydro.mass_density(*s.m),
        s.prim.pressure(*s.m),
        s.prim.specific_internal_energy(*s.m),
        s.cons.hydro.total_energy_density(*s.m),
        s.prim.sound_speed(*s.m)};

      annotation_rguard<tasks::apply_boundaries_task>(),
        sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
          *s.m,
          s.icst.bmap(*s.gt),
          scalar_v,
          std::vector{
            s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});
    }

    // Update INTEGER material field in ghost cells BEFORE flux computation
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

#endif
  // ------ ADDED NEW 1 ------- END

#ifndef IN_TEST
  // Diagnosis (only for non-test builds)
  annotation_rguard<diagnose_psi_range_1d_task>(),
    sc.execute<diagnose_psi_range_1d<D>>(flecsi::exec::on,
      *s.m,
      s.ls.psi(*s.m),
      s.prim.velocity(*s.m),
      "psi after determine_psi");
#endif

} // END of ssprk3_step function

#ifdef IN_TEST
// Explicit template instantiations for test
template void ssprk3_step<1>(flecsi::scheduler &,
  tests::level_set_utils::state<1> &);
template void ssprk3_step<2>(flecsi::scheduler &,
  tests::level_set_utils::state<2> &);
template void ssprk3_step<3>(flecsi::scheduler &,
  tests::level_set_utils::state<3> &);
#else
template<std::size_t D>
void
update_ls(control_policy<state, D> & cp) {
  ssprk3_step(cp);
}

inline control<state, 1>::action<update_ls<1>, cp::update_ls> update_ls_1d;
inline control<state, 2>::action<update_ls<2>, cp::update_ls> update_ls_2d;
inline control<state, 3>::action<update_ls<3>, cp::update_ls> update_ls_3d;
#endif // IN_TEST

} // namespace hard