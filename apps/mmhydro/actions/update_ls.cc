#include "state.hh"

#include "gfm/tasks/ls_core.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard {

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

#ifndef IN_TEST
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();
#endif

#ifdef IN_TEST
  sc.execute<tests::level_set_utils::copy_test_vel<D>>(
    flecsi::exec::on, *s.m, s.velocity(*s.m), s.t(*s.gt));
#endif // IN_TEST

  // 0) Freeze material for this entire time step: Lfrozen = L^n
  sc.execute<freeze_material<D>>(
    flecsi::exec::on, *s.m, s.ls.material(*s.m), s.ls.frozen_material(*s.m));

  // Stage 1: k = L(psi^n, L^n); psi1 = psi^n + dt*k
  sc.execute<rhs_upwind<D>>(flecsi::exec::on,
    *s.m,
    s.ls.psi(*s.m),
    s.prim.velocity(*s.m),
    s.ls.frozen_material(*s.m),
    s.ls.rhs(*s.m));

  sc.execute<advance_psi_s1<D, 1>>(flecsi::exec::on,
    *s.m,
    s.ls.psi(*s.m), // psi^n
    s.ls.rhs(*s.m), // k(psi^n)
    s.ls.psi1(*s.m), // psi1
    s.dt(*s.gt));

#ifdef IN_TEST
  sc.execute<tests::level_set_utils::copy_test_vel<D>>(
    flecsi::exec::on, *s.m, s.prim.velocity(*s.m), s.t(*s.gt));
#endif // IN_TEST

  // Stage 2: k = L(psi1, L^n); psi2 = 0.75*psi^n + 0.25*(psi1 + dt*k)
  sc.execute<rhs_upwind<D>>(flecsi::exec::on,
    *s.m,
    s.ls.psi1(*s.m),
    s.prim.velocity(*s.m),
    s.ls.frozen_material(*s.m),
    s.ls.rhs(*s.m));

  sc.execute<advance_psi<D, 2>>(flecsi::exec::on,
    *s.m,
    s.ls.psi(*s.m), // psi^n
    s.ls.rhs(*s.m), // k(psi1)
    s.ls.psi1(*s.m), // psi_prev = psi1
    s.ls.psi2(*s.m), // psi2
    s.dt(*s.gt));

#ifdef IN_TEST
  sc.execute<tests::level_set_utils::copy_test_vel<D>>(
    flecsi::exec::on, *s.m, s.prim.velocity(*s.m), s.t(*s.gt));
#endif // IN_TEST

  // Stage 3: k = L(psi2, L^n); psi3 = 1/3*psi^n + 2/3*(psi2 + dt*k)
  sc.execute<rhs_upwind<D>>(flecsi::exec::on,
    *s.m,
    s.ls.psi2(*s.m),
    s.prim.velocity(*s.m),
    s.ls.frozen_material(*s.m),
    s.ls.rhs(*s.m));

  sc.execute<advance_psi<D, 3>>(flecsi::exec::on,
    *s.m,
    s.ls.psi(*s.m), // psi^n
    s.ls.rhs(*s.m), // k(psi2)
    s.ls.psi2(*s.m), // psi_prev = psi2
    s.ls.psi3(*s.m), // psi3
    s.dt(*s.gt));

  // Algorithm 1 post-process: (psi3, L^n) -> (psi^{n+1}, L^{n+1})
  sc.execute<determine_psi<D>>(flecsi::exec::on,
    *s.m,
    s.ls.psi3(*s.m), // frozen_psi = psi3
    s.ls.psi(*s.m), // output psi^{n+1}
    s.ls.frozen_material(*s.m), // L^n
    s.ls.material(*s.m)); // output L^{n+1}
}

template<std::size_t D>
void
update_ls(control_policy<state, D> & cp) {

  ssprk3_step(cp);

} // update_ls

inline control<state, 1>::action<update_ls<1>, cp::update_ls> update_ls_1d;
inline control<state, 2>::action<update_ls<2>, cp::update_ls> update_ls_2d;
inline control<state, 3>::action<update_ls<3>, cp::update_ls> update_ls_3d;

} // namespace hard
