#ifndef HARD_MODULES_GFM_TASKS_EXTRAPOLATION_HH
#define HARD_MODULES_GFM_TASKS_EXTRAPOLATION_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/limiter.hh"
#include "../../spec/mesh.hh"
#include "../state.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::extrapolation {

// Upwind derivative with TVD blend between 1st-order and 2nd-order upwind.
// Uses minmod limiter.
// NOTE: In Aslam's paper, they are using minmod limter in
// "psi_mm(r = max(0,min(1,r))" where r is limiter ratio form.
// This should be identical to minmod on slopes for monotone data which HARD
// is already equipped with

// We follow cell-centered index:
// i-2       i-1        i        i+1       i+2
// u_{i-2}  u_{i-1}    u_i      u_{i+1}   u_{i+2}

// For sten_size == 5:
//  u[0]     u[1]      u[2]      u[3]      u[4]

// For sten_size == 3:
//           u[0]      u[1]      u[2]
//

template<std::size_t sten_size = 5>
double
gradu_1d(std::array<double, sten_size> u, double a, double dx) {
  static_assert(
    sten_size == 3 || sten_size == 5, "Stencil size must be 3 or 5");

  using limiter = spec::limiters::genminmod;

  if(a >= 0.0) {
    if constexpr(sten_size == 3) {
      // 1st order upwind
      return (u[1] - u[0]) / dx;
    }
    else {
      // 2nd-order upwind (biased)
      const double d1 = (u[2] - u[1]) / dx;
      const double d2 = (3.0 * u[2] - 4.0 * u[1] + u[0]) / (2.0 * dx);
      const double psi = limiter(u[0], u[1], u[2]);
      return d1 + psi * (d2 - d1);
    }
  }
  else {
    if constexpr(sten_size == 3) {
      // 1st order downwind
      return (u[2] - u[1]) / dx;
    }
    else {
      // 2nd-order downwind (biased)
      const double d1 = (u[2] - u[1]) / dx;
      const double d2 = (-3.0 * u[2] + 4.0 * u[3] - u[4]) / (2.0 * dx);
      const double psi = limiter(u[2], u[4], u[4]);
      return d1 + psi * (d2 - d1);
    }
  }
}

// gradients and unit normals (for GFM & extrapolation)
// centered gradient: dpsi/dx, dpsi/dy, dpsi/dz at cell centers
template<std::size_t D>
void
gradient_centered(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<rw, ro> dp_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto dp = m.template mdcolex<is::cells>(dp_a);

  const double div_dx = 1 / (2.0 * m.template delta<ax::x>());

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dp(i).x = (psi(i + 1) - psi(i - 1)) * div_dx;
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double div_dy = 1 / (2.0 * m.template delta<ax::y>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      dp(i, j).x = (psi(i + 1, j) - psi(i - 1, j)) * div_dx;
      dp(i, j).y = (psi(i, j + 1) - psi(i, j - 1)) * div_dy;
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double div_dy = 1 / (2.0 * m.template delta<ax::y>());
    const double div_dz = 1 / (2.0 * m.template delta<ax::z>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      dp(i, j, k).x = (psi(i + 1, j, k) - psi(i - 1, j, k)) * div_dx;
      dp(i, j, k).y = (psi(i, j + 1, k) - psi(i, j - 1, k)) * div_dy;
      dp(i, j, k).z = (psi(i, j, k + 1) - psi(i, j, k - 1)) * div_dz;
    };
  }
}

// Unit normal computation: n = \nabla psi / |\nabla psi|
template<std::size_t D>
void
unit_normals(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> dp_a,
  typename field<double>::template accessor<rw> normal_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto dp = m.template mdcolex<is::cells>(dp_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);

  constexpr double eps = 1e-12;

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      normal(i).x = dp(i).x / (dp(i).normal() + eps);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(dp,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      const double inv_mag = 1.0 / (dp(i, j).normal() + eps);

      normal(i, j).x = dp(i, j).x * inv_mag;
      normal(i, j).y = dp(i, j).y * inv_mag;
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(dp,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      const double inv_mag = 1.0 / (dp(i, j, k).normal() + eps);

      normal(i, j, k).x = dp(i, j, k).x * inv_mag;
      normal(i, j, k).y = dp(i, j, k).y * inv_mag;
      normal(i, j, k).z = dp(i, j, k).z * inv_mag;
    };
  }
}

// constant-normal extrapolation utility
// Extends scalar q off the interface (one forward Euler step of
// \partial tau q  + sign(\psi) n·\nabla q = 0 ).
template<std::size_t D>
void
extend_constant_normal(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<double>::template accessor<ro, ro> normal_a,
  double dtau,
  typename field<double>::template accessor<ro> q_scalar_a,
  typename field<double>::template accessor<rw> q_new_scalar_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);
  auto q = m.template mdcolex<is::cells>(q_scalar_a);
  auto q_new = m.template mdcolex<is::cells>(q_new_scalar_a);

  const double inv_dx = 1 / m.template delta<ax::x>();

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {

      // Upwind neighbor based on sign of nx
      if(normal(i).x > 0.0) {
        q_new(i) = q(i) - dtau * copysign(1.0, psi(i)) * normal(i).x *
                            (q(i) - q(i - 1)) * inv_dx;
      }
      else {
        q_new(i) = q(i);
      }
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(q,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double inv_dy = 1 / m.template delta<ax::y>();

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      // upwind neighbors (component-wise)
      const double dqx = (normal(i, j).x > 0.0) ? (q(i, j) - q(i - 1, j)) : 0.0;
      const double dqy = (normal(i, j).y > 0.0) ? (q(i, j) - q(i, j - 1)) : 0.0;

      q_new(i, j) = q(i, j) - dtau * copysign(1.0, psi(i, j)) *
                                (normal(i, j).x * dqx * inv_dx +
                                  normal(i, j).y * dqy * inv_dy);
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(q,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double inv_dy = m.template delta<ax::y>();
    const double inv_dz = m.template delta<ax::z>();

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      const double dqx =
        (normal(i, j, k).x > 0.0) ? (q(i, j, k) - q(i - 1, j, k)) : 0.0;
      const double dqy =
        (normal(i, j, k).y > 0.0) ? (q(i, j, k) - q(i, j - 1, k)) : 0.0;
      const double dqz =
        (normal(i, j, k).z > 0.0) ? (q(i, j, k) - q(i, j, k - 1)) : 0.0;

      q_new(i, j, k) = q(i, j, k) - dtau * copysign(1.0, psi(i, j, k)) *
                                      (normal(i, j, k).x * dqx * inv_dx +
                                        normal(i, j, k).y * dqy * inv_dy +
                                        normal(i, j, k).z * dqz * inv_dz);
      ;
    };
  }
}

// ----------------------------
// RHS: R(u) = - H(psi) * (n · grad(u))  in psi>0, else 0.
// grad(u) computed with 2nd-order upwind + MinMod for each component.
// ----------------------------

template<std::size_t D>
void
rhs_constant_extrapolation_minmod(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> q_scalar_a,
  typename field<double>::template accessor<ro> psi_a,
  typename field<double>::template accessor<ro> normal_a,
  typename field<double>::template accessor<rw> rhs_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto q = m.template mdcolex<is::cells>(q_scalar_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);

  const double dx = m.template delta<ax::x>();

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {

      if(psi(i) <= 0.0) {
        rhs(i) = 0.0;
        return;
      }

      const double ax = normal(i).x;

      rhs(i) =
        -ax * gradu_1d({q(i - 2), q(i - 1), q(i), q(i + 1), q(i + 2)}, ax, dx);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(q,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double dy = m.template delta<ax::y>();

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      if(psi(i, j) <= 0.0) {
        rhs(i, j) = 0.0;
        return;
      }

      const double ax = normal(i, j).x;
      const double ay = normal(i, j).y;

      // clang-format off
      rhs(i, j) =
        -(ax * gradu_1d(
          {q(i - 2, j), q(i - 1, j), q(i, j), q(i + 1, j), q(i + 2, j)},
          ax, dx) +
        ay * gradu_1d(
          {q(i, j - 2), q(i, j - 1), q(i, j), q(i, j + 1), q(i, j + 2)},
          ay, dy));
      // clang-format on
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(q,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double dy = m.template delta<ax::y>();
    const double dz = m.template delta<ax::z>();

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      if(psi(i, j, k) <= 0.0) {
        rhs(i, j, k) = 0.0;
        return;
      }

      const double ax = normal(i, j, k).x;
      const double ay = normal(i, j, k).y;
      const double az = normal(i, k, k).z;

      // clang-format off
      rhs(i, j, k) =
        -(ax * gradu_1d(
          {q(i - 2, j, k),
            q(i - 1, j, k),
            q(i, j, k),
            q(i + 1, j, k),
            q(i + 2, j, k)},
          ax, dx) +
        ay * gradu_1d(
          {q(i, j - 2, k),
            q(i, j - 1, k),
            q(i, j, k),
            q(i, j + 1, k),
            q(i, j + 2, k)},
          ay, dy) +
        az * gradu_1d(
          {q(i, j, k - 2),
            q(i, j, k - 1),
            q(i, j, k),
            q(i, j, k + 1),
            q(i, j, k + 2)},
          az, dz));
      // clang-format on
    };
  }
}

// RK2(Heun) advance for constant extrapolation
// STEP=1: q_stage = q + dt * rhs
// STEP=2: q_out   = 0.5*( q + q_stage + dt*rhs_stage )
template<std::size_t D, std::size_t STEP>
void
advance_q_rk2(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> psi_a,
  typename field<double>::template accessor<ro> q_scalar_known_a,
  typename field<double>::template accessor<ro> q_scalar_a, // base q^n
  typename field<double>::template accessor<ro> rhs_a, // rhs at base or stage
  typename field<double>::template accessor<ro>
    q_scalar_stage_a, // q^(1) (only used in STEP=2)
  typename field<double>::template accessor<wo>
    q_scalar_new_a, // q^(1) for STEP=1, q^(n+1) for STEP=2
  const double dt) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  static_assert(STEP == 1 || STEP == 2, "STEP is 1 or 2");

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto q_known = m.template mdcolex<is::cells>(q_scalar_known_a);
  auto q = m.template mdcolex<is::cells>(q_scalar_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);
  auto q_stage = m.template mdcolex<is::cells>(q_scalar_stage_a);
  auto q_new = m.template mdcolex<is::cells>(q_scalar_new_a);

  auto is_known = [](double ph) { return (ph <= 0.0); };

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(is_known(psi(i))) {
        q_new(i) = q_known(i);
        return;
      }

      if constexpr(STEP == 1) {
        q_new(i) = q(i) + dt * rhs(i); // q^(1)
      }
      else {
        // rhs is rhs(q_stage) here
        q_new(i) = 0.5 * (q(i) + q_stage(i) + dt * rhs(i));
      }
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(q,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      if(is_known(psi(i, j))) {
        q_new(i, j) = q_known(i, j);
        return;
      }

      if constexpr(STEP == 1) {
        q_new(i, j) = q(i, j) + dt * rhs(i, j);
      }
      else {
        q_new(i, j) = 0.5 * (q(i, j) + q_stage(i, j) + dt * rhs(i, j));
      }
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(q,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      if(is_known(psi(i, j, k))) {
        q_new(i, j, k) = q_known(i, j, k);
        return;
      }

      if constexpr(STEP == 1) {
        q_new(i, j, k) = q(i, j, k) + dt * rhs(i, j, k);
      }
      else {
        q_new(i, j, k) =
          0.5 * (q(i, j, k) + q_stage(i, j, k) + dt * rhs(i, j, k));
      }
    };
  }
}

} // namespace hard::tasks::extrapolation

namespace hard::extrapolation {

template<std::size_t D>
void
rk2_constant_extrapolation_step(control_policy<hard::state, D> & cp) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  // Stage 1: rhs(q^n) -> rhs
  sc.execute<tasks::extrapolation::rhs_constant_extrapolation_minmod<D>>(
    flecsi::exec::on,
    *s.m,
    s.extp.q_scalar(*s.m),
    s.ls.psi(*s.m),
    s.extp.normal(*s.m),
    s.extp.q_rhs(*s.m));

  // q1 = q^n + dt * rhs(q^n)
  sc.execute<tasks::extrapolation::advance_q_rk2<D, 1>>(flecsi::exec::on,
    *s.m,
    s.ls.psi(*s.m),
    s.extp.q_scalar_known(*s.m),
    s.extp.q_scalar(*s.m),
    s.extp.q_rhs(*s.m),
    s.extp.q_scalar2(*s.m), // dummy (unused in STEP=1)
    s.extp.q_sclar1(*s.m), // output stage-1
    s.dt(s.gt)); // dtau in the origianl work

  // Stage 2: rhs(q1) -> rhs1
  sc.execute<tasks::extrapolation::rhs_constant_extrapolation_minmod<D>>(
    flecsi::exec::on,
    *s.m,
    s.extp.q_scalar1(*s.m),
    s.ls.psi(*s.m),
    s.extp.normal(*s.m),
    s.extp.q_rhs1(*s.m));

  // q^{n+1} = 0.5*( q^n + q1 + dt*rhs(q1) )
  sc.execute<tasks::extrapolation::advance_q_rk2<D, 2>>(flecsi::exec::on,
    *s.m,
    s.ls.psi(*s.m),
    s.extp.q_scalar_known(*s.m),
    s.extp.q_scalar(*s.m), // q^n
    s.extp.q_rhs1(*s.m), // rhs(q1)
    s.extp.q_scalar1(*s.m), // q_stage
    s.extp.q_scalar2(*s.m), // q_new
    s.dt(s.gt)); // dtau in original work

  // Swap
  // sc.execute<hard::tasks::util::copy_scal_field<dm::quantities, D, double>>(
  //   flecsi::exec::on, *s.m, s.expt.q_scalar2(*s.m), s.expt.q_scalar(*s.m));
}

} // namespace hard::extrapolation

#endif
