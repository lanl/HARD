#ifndef HARD_MODULES_GFM_TASKS_EXTRAPOLATION_HH
#define HARD_MODULES_GFM_TASKS_EXTRAPOLATION_HH

#include "../../../apps/mmhydro/options.hh"
#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/limiter.hh"
#include "../../spec/mesh.hh"
#include "../../spec/tasks/io.hh"
#include "../state.hh"
#include "utils.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace hard::tasks::extrapolation {

// Default pseudo-time controls for the Aslam (JCP 2003) constant
// extrapolation routine used to fill ghost states.
inline constexpr double extp_cfl = 1.0; // we have a linear PDE problem
inline constexpr double extp_residual_tol =
  1.0; // bound on max SciPy-scaled update in band

template<std::size_t D, typename CP>
inline void
debug_output_extrapolation_t0(CP & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  spec::config_py config(opt::config.value());

  const bool debug_output = config.contains("debug_extrapolation_t0") &&
                            config["debug_extrapolation_t0"].cast<bool>();

  if(!debug_output) {
    return;
  }

  static bool wrote_debug_output = false;
  if(wrote_debug_output) {
    return;
  }

  flog(info) << "Writing debug extrapolation output at initial stage.\n";

  sc.execute<tasks::util::material_to_scalar_output<D>>(
    flecsi::exec::on, *s.m, s.ls.material(*s.m), s.ls.material_out(*s.m));

  sc.execute<tasks::util::extp_rho_to_scalar_output<D>>(flecsi::exec::on,
    *s.m,
    s.extp.q_ext_plus(*s.m),
    s.extp.q_ext_plus_rho_out(*s.m));

  sc.execute<tasks::util::extp_rho_to_scalar_output<D>>(flecsi::exec::on,
    *s.m,
    s.extp.q_ext_minus(*s.m),
    s.extp.q_ext_minus_rho_out(*s.m));

  sc.execute<tasks::util::vec_components_to_scalar_output<D>>(flecsi::exec::on,
    *s.m,
    s.extp.dp(*s.m),
    s.extp.dp_x_out(*s.m),
    s.extp.dp_y_out(*s.m),
    s.extp.dp_z_out(*s.m),
    s.extp.dp_mag_out(*s.m));

  sc.execute<tasks::util::vec_components_to_scalar_output<D>>(flecsi::exec::on,
    *s.m,
    s.extp.normal(*s.m),
    s.extp.normal_x_out(*s.m),
    s.extp.normal_y_out(*s.m),
    s.extp.normal_z_out(*s.m),
    s.extp.normal_mag_out(*s.m));

  auto lm = flecsi::data::launch::make(sc, *s.m);

  auto scalar_v = std::vector{std::make_tuple(s.cons.hydro.mass_density(lm),
                                std::string{"density_before_gfm"}),
    std::make_tuple(
      s.extp.q_ext_plus_rho_out(lm), std::string{"q_ext_plus_density"}),
    std::make_tuple(
      s.extp.q_ext_minus_rho_out(lm), std::string{"q_ext_minus_density"}),
    std::make_tuple(s.extp.extp_psi(lm), std::string{"extp_psi"}),
    std::make_tuple(s.ls.psi(lm), std::string{"ls_psi"}),
    std::make_tuple(s.ls.material_out(lm), std::string{"material_id"}),
    std::make_tuple(s.extp.dp_x_out(lm), std::string{"extp_dp_x"}),
    std::make_tuple(s.extp.dp_y_out(lm), std::string{"extp_dp_y"}),
    std::make_tuple(s.extp.dp_z_out(lm), std::string{"extp_dp_z"}),
    std::make_tuple(s.extp.dp_mag_out(lm), std::string{"extp_dp_mag"}),
    std::make_tuple(s.extp.normal_x_out(lm), std::string{"extp_normal_x"}),
    std::make_tuple(s.extp.normal_y_out(lm), std::string{"extp_normal_y"}),
    std::make_tuple(s.extp.normal_z_out(lm), std::string{"extp_normal_z"}),
    std::make_tuple(s.extp.normal_mag_out(lm), std::string{"extp_normal_mag"})};

  if(s.output_method == spec::om::vti) {
    flecsi::execute<tasks::io::vtk<D>, flecsi::mpi>(flecsi::exec::on,
      spec::io::name{"debug_extrapolation_t0"},
      s.t(*s.gt),
      lm,
      scalar_v,
      std::vector{
        std::make_tuple(s.prim.velocity(lm), std::string{"velocity"})});
  }
#ifdef HARD_ENABLE_HDF5
  else if(s.output_method == spec::om::xdmf) {
    flecsi::execute<tasks::io::xdmf<D>, flecsi::mpi>(flecsi::exec::on,
      spec::io::name{"debug_extrapolation_t0"},
      s.t(*s.gt),
      lm,
      scalar_v,
      std::vector{
        std::make_tuple(s.prim.velocity(lm), std::string{"velocity"})});
  }
#endif
  else {
    flecsi::execute<tasks::io::csv<D>, flecsi::mpi>(flecsi::exec::on,
      spec::io::name{"debug_extrapolation_t0"},
      s.t(*s.gt),
      lm,
      scalar_v,
      std::vector{
        std::make_tuple(s.prim.velocity(lm), std::string{"velocity"})});
  }

  wrote_debug_output = true;
}

FLECSI_INLINE_TARGET bool
same_extp_side(double a, double b) noexcept {
  return (a <= 0.0 && b <= 0.0) || (a > 0.0 && b > 0.0);
}

template<typename Fetch>
FLECSI_INLINE_TARGET double
side_aware_derivative_1d(Fetch fetch, double pi, double dx) noexcept {
  // Always need immediate neighbors to decide what stencil is possible.
  const double pim1 = fetch(-1);
  const double pip1 = fetch(1);
  const bool left1_ok = same_extp_side(pi, pim1);
  const bool right1_ok = same_extp_side(pi, pip1);
  // Best case: centered second-order derivative.
  // Only two neighbor reads.
  if(left1_ok && right1_ok) {
    return (pip1 - pim1) / (2.0 * dx);
  }
  // Left-biased derivative.
  // Only now read i-2.
  if(left1_ok) {
    const double pim2 = fetch(-2);
    const bool left2_ok = same_extp_side(pi, pim2);
    if(left2_ok) {
      return (3.0 * pi - 4.0 * pim1 + pim2) / (2.0 * dx);
    }
    return (pi - pim1) / dx;
  }
  // Right-biased derivative.
  // Only now read i+2.
  if(right1_ok) {
    const double pip2 = fetch(2);
    const bool right2_ok = same_extp_side(pi, pip2);
    if(right2_ok) {
      return (-3.0 * pi + 4.0 * pip1 - pip2) / (2.0 * dx);
    }
    return (pip1 - pi) / dx;
  }
  // No same-side neighbor available.
  return 0.0;
}

struct linear_axis_stencil {
  double w[2];
  int idx[2];
  int n;
};

FLECSI_INLINE_TARGET linear_axis_stencil
make_linear_axis_stencil(double x) noexcept {
  constexpr double eps = 1.0e-14;
  linear_axis_stencil s{};
  const int i0 = static_cast<int>(std::floor(x));
  const double t = x - static_cast<double>(i0);
  // x is effectively exactly i0.
  if(t <= eps) {
    s.idx[0] = i0;
    s.w[0] = 1.0;
    s.n = 1;
    return s;
  }
  // x is effectively exactly i0 + 1.
  if((1.0 - t) <= eps) {
    s.idx[0] = i0 + 1;
    s.w[0] = 1.0;
    s.n = 1;
    return s;
  }
  // General case: two nonzero linear weights.
  s.idx[0] = i0;
  s.idx[1] = i0 + 1;
  s.w[0] = 1.0 - t;
  s.w[1] = t;
  s.n = 2;
  return s;
}

template<typename ScalarAccessor>
FLECSI_INLINE_TARGET double
sample_bilinear_nonzero(const ScalarAccessor & q,
  double xi,
  double yj) noexcept {
  const auto xs = make_linear_axis_stencil(xi);
  const auto ys = make_linear_axis_stencil(yj);
  double value = 0.0;
#pragma unroll 2
  for(int b = 0; b < ys.n; ++b) {
    for(int a = 0; a < xs.n; ++a) {
      const double w = xs.w[a] * ys.w[b];
      // w is nonzero by construction.
      value += w * q(xs.idx[a], ys.idx[b]);
    }
  }
  return value;
}

template<typename ScalarAccessor>
FLECSI_INLINE_TARGET double
sample_trilinear_nonzero(const ScalarAccessor & q,
  double xi,
  double yj,
  double zk) noexcept {
  const auto xs = make_linear_axis_stencil(xi);
  const auto ys = make_linear_axis_stencil(yj);
  const auto zs = make_linear_axis_stencil(zk);
  double value = 0.0;
#pragma unroll 2
  for(int c = 0; c < zs.n; ++c) {
#pragma unroll 2
    for(int b = 0; b < ys.n; ++b) {
      for(int a = 0; a < xs.n; ++a) {
        const double w = xs.w[a] * ys.w[b] * zs.w[c];
        // w is nonzero by construction.
        value += w * q(xs.idx[a], ys.idx[b], zs.idx[c]);
      }
    }
  }
  return value;
}

template<std::size_t D, typename QAccessor>
FLECSI_INLINE_TARGET double
directional_derivative_normal_1st(const QAccessor & q,
  std::size_t qid,
  const vec<D> & n,
  int i,
  int j,
  int k,
  double dx,
  double dy = 1.0,
  double dz = 1.0) noexcept {
  constexpr double eps = 1.0e-14;
  const double ax = n.x() / dx;
  double inv_ds = std::abs(ax);
  double ay = 0.0;
  double az = 0.0;
  if constexpr(D >= 2) {
    ay = n.y() / dy;
    inv_ds += std::abs(ay);
  }
  if constexpr(D >= 3) {
    az = n.z() / dz;
    inv_ds += std::abs(az);
  }
  if(inv_ds <= eps) {
    return 0.0;
  }

  const double ds = 1.0 / inv_ds;
  const double xi = static_cast<double>(i) - ds * ax;

  if constexpr(D == 1) {
    const double q0 = q(i)[qid];
    const auto xs = make_linear_axis_stencil(xi);
    double qm1 = 0.0;
    for(int a = 0; a < xs.n; ++a) {
      const int ii = xs.idx[a];
      const double val = (ii == i) ? q0 : q(ii)[qid];
      qm1 += xs.w[a] * val;
    }
    return (q0 - qm1) * inv_ds;
  }
  else if constexpr(D == 2) {
    const double q0 = q(i, j)[qid];
    const double yj = static_cast<double>(j) - ds * ay;
    const auto fetch = [&](int ii, int jj) {
      return (ii == i && jj == j) ? q0 : q(ii, jj)[qid];
    };
    const double qm1 = sample_bilinear_nonzero(fetch, xi, yj);
    return (q0 - qm1) * inv_ds;
  }
  else {
    const double q0 = q(i, j, k)[qid];
    const double yj = static_cast<double>(j) - ds * ay;
    const double zk = static_cast<double>(k) - ds * az;
    const auto fetch = [&](int ii, int jj, int kk) {
      return (ii == i && jj == j && kk == k) ? q0 : q(ii, jj, kk)[qid];
    };
    const double qm1 = sample_trilinear_nonzero(fetch, xi, yj, zk);
    return (q0 - qm1) * inv_ds;
  }
}

struct gradient_centered_task : annotation_task_region {
  inline static const std::string name{"gradient_centered"};
};

// gradients and unit normals (for GFM & extrapolation)
// centered gradient: dpsi/dx, dpsi/dy, dpsi/dz at cell centers
template<std::size_t D>
void
gradient_centered(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<wo, na> dp_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto dp = m.template mdcolex<is::cells>(dp_a);
  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double psi_c = psi(i);
      dp(i).x() =
        side_aware_derivative_1d([&](int off) { return psi(i + off); },
          psi_c,
          m.template delta<ax::x>());
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      const double psi_c = psi(i, j);
      dp(i, j).x() =
        side_aware_derivative_1d([&](int off) { return psi(i + off, j); },
          psi_c,
          m.template delta<ax::x>());
      dp(i, j).y() =
        side_aware_derivative_1d([&](int off) { return psi(i, j + off); },
          psi_c,
          m.template delta<ax::y>());
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      const double psi_c = psi(i, j, k);
      dp(i, j, k).x() =
        side_aware_derivative_1d([&](int off) { return psi(i + off, j, k); },
          psi_c,
          m.template delta<ax::x>());
      dp(i, j, k).y() =
        side_aware_derivative_1d([&](int off) { return psi(i, j + off, k); },
          psi_c,
          m.template delta<ax::y>());
      dp(i, j, k).z() =
        side_aware_derivative_1d([&](int off) { return psi(i, j, k + off); },
          psi_c,
          m.template delta<ax::z>());
    };
  }
}
struct unit_cell_normals_task : annotation_task_region {
  inline static const std::string name{"unit_cell_normals"};
};
// Unit normal computation: n = \nabla psi / |\nabla psi|
template<std::size_t D>
void
unit_cell_normals(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<vec<D>>::template accessor<ro, ro> dp_a,
  typename field<vec<D>>::template accessor<wo, na> normal_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto dp = m.template mdcolex<is::cells>(dp_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);

  constexpr double eps = 1e-12;

  auto normalize = [](const vec<D> & g) {
    vec<D> n(0.0);
    const double mag = g.norm();
    if(mag > eps) {
      n = g / mag;
    }
    return n;
  };

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      normal(i) = dp(i) / (dp(i).norm() + eps);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(dp,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      normal(i, j) = dp(i, j) / (dp(i, j).norm() + eps);
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(dp,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      normal(i, j, k) = dp(i, j, k) / (dp(i, j, k).norm() + eps);
    };
  }
}
// -------------------------------------------------------------------
// RHS: R(u) = - H(psi) * (n · grad(u))  in psi>0, else 0.
// grad(u) computed with 1st-order upwind for each component.
// -------------------------------------------------------------------

struct rhs_constant_extrapolation_task : annotation_task_region {
  inline static const std::string name{"rhs_constant_extrapolation"};
};

template<std::size_t D>
void
rhs_constant_extrapolation(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro>
    q_scalar_a,
  typename field<double>::template accessor<ro, na> psi_a,
  typename field<vec<D>>::template accessor<ro, na> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    rhs_a) noexcept {

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
      }
      else {
        const auto n = normal(i);
        for(std::size_t qid = 0; qid < q(i).size(); ++qid) {
          rhs(i)[qid] =
            -directional_derivative_normal_1st<D>(q, qid, n, i, 0, 0, dx);
        }
      }
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
      }
      else {
        const auto n = normal(i, j);
        for(std::size_t qid = 0; qid < q(i, j).size(); ++qid) {
          rhs(i, j)[qid] =
            -directional_derivative_normal_1st<D>(q, qid, n, i, j, 0, dx, dy);
        }
      }
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
      }
      else {
        const auto n = normal(i, j, k);
        for(std::size_t qid = 0; qid < q(i, j, k).size(); ++qid) {
          rhs(i, j, k)[qid] = -directional_derivative_normal_1st<D>(
            q, qid, n, i, j, k, dx, dy, dz);
        }
      }
    };
  }
}

struct extp_interface_convergence_update_task : annotation_task_region {
  inline static const std::string name{"extp_interface_convergence_update"};
};

struct extp_interface_convergence_task : annotation_task_region {
  inline static const std::string name{"extp_interface_convergence"};
};

template<std::size_t D>
double
extp_interface_convergence_update(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_old_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_new_a,
  typename field<double>::template accessor<ro, ro> psi_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;
  namespace fold = flecsi::exec::fold;

  auto q_old = m.template mdcolex<is::cells>(q_old_a);
  auto q_new = m.template mdcolex<is::cells>(q_new_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);

  // Number of physical packed quantities used in the convergence check.
  // Exclude the trailing metadata slot n_comps().
  constexpr std::size_t N_Q = hard::tasks::util::extp_q<D>::N_q - 1;

  // Tuning constants:
  //   rtol       : relative tolerance
  //   atol_frac  : absolute tol based on the max of each components of q
  //   tiny       : guards against all-zero fields (q_ref == 0)
  constexpr double rtol = 1e-2;
  constexpr double atol_frac = 1e-3;
  constexpr double tiny = 1e-16;

  // SciPy-style scaled update check in the interface band:
  //   scale_qid = atol[qid] + rtol * |q^{n+1}_qid|
  //   e_qid     = |q^{n+1}_qid - q^n_qid| / scale_qid
  //   e_cell    = max_qid(e_qid)
  // atol[qid] is derived from the band-wide max of |q^{n+1}_qid|.
  std::array<double, N_Q> q_ref{};
  std::array<double, N_Q> atol{};

  if constexpr(D == 1) {
    auto in_band = [=](auto i) {
      return psi(i) > 0.0 && (psi(i - 1) <= 0.0 || psi(i + 1) <= 0.0);
    };

    // --- Pass 1: per-variable band-max reference scale ---
    for(std::size_t qid = 0; qid < N_Q; ++qid) {
      q_ref[qid] = s.executor().reduceall(
        i, up, (m.template cells<ax::x, dm::quantities>()), fold::max, double) {
        up(in_band(i) ? std::abs(q_new(i)[qid]) : 0.0);
      };
    }
    for(std::size_t qid = 0; qid < N_Q; ++qid)
      atol[qid] = std::max(atol_frac * q_ref[qid], tiny);

    auto scaled_update = [=](const auto & qo, const auto & qn) {
      double out = 0.0;
      for(std::size_t qid = 0; qid < N_Q; ++qid) {
        const double diff = std::abs(qn[qid] - qo[qid]);
        const double denom = atol[qid] + rtol * std::abs(qn[qid]);
        out = std::max(out, diff / denom);
      }
      return out;
    };

    // --- Pass 2: scaled-update reduction ---
    return s.executor().reduceall(
      i, up, (m.template cells<ax::x, dm::quantities>()), fold::max, double) {
      up(in_band(i) ? scaled_update(q_old(i), q_new(i)) : 0.0);
    };
  }
  else if constexpr(D == 2) {
    auto in_band = [=](auto i, auto j) {
      return psi(i, j) > 0.0 && (psi(i - 1, j) <= 0.0 || psi(i + 1, j) <= 0.0 ||
                                  psi(i, j - 1) <= 0.0 || psi(i, j + 1) <= 0.0);
    };

    auto mdpolicy_qq = get_mdiota_policy(q_new,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    // --- Pass 1: per-variable band-max reference scale ---
    for(std::size_t qid = 0; qid < N_Q; ++qid) {
      q_ref[qid] =
        s.executor().reduceall(ji, up, mdpolicy_qq, fold::max, double) {
        auto [j, i] = ji;
        up(in_band(i, j) ? std::abs(q_new(i, j)[qid]) : 0.0);
      };
    }
    for(std::size_t qid = 0; qid < N_Q; ++qid)
      atol[qid] = std::max(atol_frac * q_ref[qid], tiny);

    auto scaled_update = [=](const auto & qo, const auto & qn) {
      double out = 0.0;
      for(std::size_t qid = 0; qid < N_Q; ++qid) {
        const double diff = std::abs(qn[qid] - qo[qid]);
        const double denom = atol[qid] + rtol * std::abs(qn[qid]);
        out = std::max(out, diff / denom);
      }
      return out;
    };

    // --- Pass 2: scaled-update reduction ---
    return s.executor().reduceall(ji, up, mdpolicy_qq, fold::max, double) {
      auto [j, i] = ji;
      up(in_band(i, j) ? scaled_update(q_old(i, j), q_new(i, j)) : 0.0);
    };
  }
  else {
    auto in_band = [=](auto i, auto j, auto k) {
      return psi(i, j, k) > 0.0 &&
             (psi(i - 1, j, k) <= 0.0 || psi(i + 1, j, k) <= 0.0 ||
               psi(i, j - 1, k) <= 0.0 || psi(i, j + 1, k) <= 0.0 ||
               psi(i, j, k - 1) <= 0.0 || psi(i, j, k + 1) <= 0.0);
    };

    auto mdpolicy_qqq = get_mdiota_policy(q_new,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    // --- Pass 1: per-variable band-max reference scale ---
    for(std::size_t qid = 0; qid < N_Q; ++qid) {
      q_ref[qid] =
        s.executor().reduceall(kji, up, mdpolicy_qqq, fold::max, double) {
        auto [k, j, i] = kji;
        up(in_band(i, j, k) ? std::abs(q_new(i, j, k)[qid]) : 0.0);
      };
    }
    for(std::size_t qid = 0; qid < N_Q; ++qid)
      atol[qid] = std::max(atol_frac * q_ref[qid], tiny);

    auto scaled_update = [=](const auto & qo, const auto & qn) {
      double out = 0.0;
      for(std::size_t qid = 0; qid < N_Q; ++qid) {
        const double diff = std::abs(qn[qid] - qo[qid]);
        const double denom = atol[qid] + rtol * std::abs(qn[qid]);
        out = std::max(out, diff / denom);
      }
      return out;
    };

    // --- Pass 2: scaled-update reduction ---
    return s.executor().reduceall(kji, up, mdpolicy_qqq, fold::max, double) {
      auto [k, j, i] = kji;
      up(
        in_band(i, j, k) ? scaled_update(q_old(i, j, k), q_new(i, j, k)) : 0.0);
    };
  }
}

template<std::size_t D>
void
extp_interface_convergence(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro>,
  flecsi::future<double> max_scaled_update_f,
  int source_material,
  int steps,
  double tol) noexcept {

  const double max_scaled_update = max_scaled_update_f.get();

  if(max_scaled_update > 10 * tol) {
    flog(warn) << "Pseudo-time PDE did not reach steady state "
                  "within the interface band.\n"
               << "  - Max_scaled_update = " << max_scaled_update
               << " , while set tol = " << tol << "\n"
               << "  - Increase 'extp_steps' (currently " << steps
               << "). For 1D problems ~ 20+, 2D ~ 50+, 3D ~ 100+.\n";
  }
  else if(max_scaled_update > tol) {
    flog(info)
      << "Pseudo-time PDE partially converged, increase 'extp_steps' (now = "
      << steps << ") by a few steps.\n";
  }
}

template<std::size_t D>
FLECSI_INLINE_TARGET hard::tasks::util::extp_q<D>
q_stage1(const hard::tasks::util::extp_q<D> & q0,
  const hard::tasks::util::extp_q<D> & rhs,
  double dtau) noexcept {

  hard::tasks::util::extp_q<D> out{};

  for(std::size_t val_i = 0; val_i < out.size(); val_i++) {
    out[val_i] = q0[val_i] + dtau * rhs[val_i];
  }
  return out;
}

template<std::size_t D>
FLECSI_INLINE_TARGET hard::tasks::util::extp_q<D>
q_stage2(const hard::tasks::util::extp_q<D> & q0,
  const hard::tasks::util::extp_q<D> & q1,
  const hard::tasks::util::extp_q<D> & rhs1,
  double dtau) noexcept {

  hard::tasks::util::extp_q<D> out{};

  for(std::size_t val_i = 0; val_i < out.size(); val_i++) {
    out[val_i] = 0.5 * (q0[val_i] + q1[val_i] + dtau * rhs1[val_i]);
  }

  return out;
}

struct advance_q_rk2_s1_cfl_task : annotation_task_region {
  inline static const std::string name{"advance_q_rk2_s1_cfl"};
};

template<std::size_t D>
void
advance_q_rk2_s1_cfl(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> extp_psi_a,
  typename field<vec<D>>::template accessor<ro, na> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> q0_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> rhs_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na> q1_a,
  double extp_cfl) noexcept {

  const double dx = m.template delta<ax::x>();
  double dy = 1.0;
  double dz = 1.0;
  if constexpr(D >= 2)
    dy = m.template delta<ax::y>();
  if constexpr(D >= 3)
    dz = m.template delta<ax::z>();
  constexpr double eps = 1e-14;

  s.executor().forall(i, flecsi::util::iota_view({}, q_known_a.span().size())) {
    if(extp_psi_a(i) <= 0.0) {
      q1_a(i) = q_known_a(i);
    }
    else {
      const auto n = normal_a(i);
      double denom = std::abs(n.x()) / dx;
      if constexpr(D >= 2)
        denom += std::abs(n.y()) / dy;
      if constexpr(D >= 3)
        denom += std::abs(n.z()) / dz;

      const double dtau = extp_cfl / std::max(denom, eps);
      q1_a(i) = q_stage1<D>(q0_a(i), rhs_a(i), dtau);
    }
  };
}

struct advance_q_rk2_s2_cfl_task : annotation_task_region {
  inline static const std::string name{"advance_q_rk2_s2_cfl"};
};

// Evolving solution q^n, updated IN PLACE to q^{n+1}
template<std::size_t D>
void
advance_q_rk2_s2_cfl(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> extp_psi_a,
  typename field<vec<D>>::template accessor<ro, na> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<rw, na> q0_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    rhs1_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> q1_a,
  double extp_cfl) noexcept {

  const double dx = m.template delta<ax::x>();
  double dy = 1.0;
  double dz = 1.0;
  if constexpr(D >= 2)
    dy = m.template delta<ax::y>();
  if constexpr(D >= 3)
    dz = m.template delta<ax::z>();
  constexpr double eps = 1e-14;

  s.executor().forall(i, flecsi::util::iota_view({}, q_known_a.span().size())) {
    if(extp_psi_a(i) <= 0.0) {
      q0_a(i) = q_known_a(i);
    }
    else {
      const auto n = normal_a(i);
      double denom = std::abs(n.x()) / dx;
      if constexpr(D >= 2)
        denom += std::abs(n.y()) / dy;
      if constexpr(D >= 3)
        denom += std::abs(n.z()) / dz;

      const double dtau = extp_cfl / std::max(denom, eps);
      q0_a(i) = q_stage2<D>(q0_a(i), q1_a(i), rhs1_a(i), dtau);
    }
  };
}

struct copy_extp_q_task : annotation_task_region {
  inline static const std::string name{"copy_extp_q"};
};

template<std::size_t D>
void
copy_extp_q(flecsi::exec::accelerator s,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> src_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    dst_a) noexcept {

  s.executor().forall(i, flecsi::util::iota_view({}, src_a.span().size())) {
    dst_a(i) = src_a(i);
  };
}

struct prepare_extp_for_material_task : annotation_task_region {
  inline static const std::string name{"prepare_extp_for_material"};
};

// Prepare one material's state for constant extrapolation.
//
// source_material:
//   material ID whose primitive state is treated as the known region.
//
// extp_psi:
//   signed extrapolation level set built from source-material membership:
//     extp_psi <= 0 in source_material
//     extp_psi >  0 outside source_material
//
//   This is done without assuming a fixed global relation between material ID
//   and the sign of ls_psi. The source material is always mapped to the known
//   side by using +/- abs(ls_psi).
//
// q/q_known:
//   packed current primitive state. q is the extrapolation work field;
//   q_known is used by RK to keep the known region fixed.
template<std::size_t D>
void
prepare_extp_for_material(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, na> ls_psi_a,
  typename field<int>::template accessor<ro, na> material_a,
  typename field<double>::template accessor<ro, na> mass_density_a,
  typename field<vec<D>>::template accessor<ro, na> momentum_density_a,
  typename field<double>::template accessor<ro, na> total_energy_density_a,
  typename field<double>::template accessor<wo, na> extp_psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na> q_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_known_a,
  int source_material) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto ls_psi = m.template mdcolex<is::cells>(ls_psi_a);
  auto material = m.template mdcolex<is::cells>(material_a);
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto extp_psi = m.template mdcolex<is::cells>(extp_psi_a);
  auto q = m.template mdcolex<is::cells>(q_a);
  auto q_known = m.template mdcolex<is::cells>(q_known_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::all>())) {
      extp_psi(i) = (material(i) == source_material) ? -std::abs(ls_psi(i))
                                                     : std::abs(ls_psi(i));

      q(i).rho() = mass_density(i);
      q(i).ru() = momentum_density(i);
      q(i).re() = total_energy_density(i);
      q_known(i) = q(i);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::all>(),
      m.template cells<ax::x, dm::all>());

    s.executor().forall(ji, mdpolicy_yx) {
      auto [j, i] = ji;
      extp_psi(i, j) = (material(i, j) == source_material)
                         ? -std::abs(ls_psi(i, j))
                         : std::abs(ls_psi(i, j));

      q(i, j).rho() = mass_density(i, j);
      q(i, j).ru() = momentum_density(i, j);
      q(i, j).re() = total_energy_density(i, j);

      q_known(i, j) = q(i, j);
    };
  }
  else {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::all>(),
      m.template cells<ax::y, dm::all>(),
      m.template cells<ax::x, dm::all>());

    s.executor().forall(kji, mdpolicy_zyx) {
      auto [k, j, i] = kji;
      extp_psi(i, j, k) = (material(i, j, k) == source_material)
                            ? -std::abs(ls_psi(i, j, k))
                            : std::abs(ls_psi(i, j, k));

      q(i, j, k).rho() = mass_density(i, j, k);
      q(i, j, k).ru() = momentum_density(i, j, k);
      q(i, j, k).re() = total_energy_density(i, j, k);

      q_known(i, j, k) = q(i, j, k);
    };
  }
}

} // namespace hard::tasks::extrapolation

#endif
