#ifndef HARD_MODULES_GFM_TASKS_EXTRAPOLATION_HH
#define HARD_MODULES_GFM_TASKS_EXTRAPOLATION_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/limiter.hh"
#include "../../spec/mesh.hh"
#include "../state.hh"
#include "utils.hh"
#include <array>
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

FLECSI_INLINE_TARGET double
limiter_minmod(double r) {
  // Original minmod in Aslam's paper limiter in ratio form: max(0, min(1, r))
  if(r <= 0.0)
    return 0.0;
  return (r < 1.0) ? r : 1.0;
}

template<std::size_t sten_size = 5>
FLECSI_INLINE_TARGET double
gradu_1d(const std::array<double, sten_size> & u, double a, double dx) {
  static_assert(
    sten_size == 3 || sten_size == 5, "Stencil size must be 3 or 5");

  // TODO : Check limiter in HARD is functional as expected
  // using limiter = spec::limiters::genminmod;

  if constexpr(sten_size == 3) {
    // Plain 1st order upwind using u = {u_{i-1}, u_i, u_{i+1}}
    if(a >= 0.0)
      return (u[1] - u[0]) / dx;
    else
      return (u[2] - u[1]) / dx;
  }
  else {
    // 5-point version: u = {u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2}}
    constexpr double eps = 1e-14;

    if(a >= 0.0) {
      // 1st order upwind derivative at i
      const double d1 = (u[2] - u[1]) / dx;

      // 2nd order upwind (biased left) derivative at i
      const double d2 = (3.0 * u[2] - 4.0 * u[1] + u[0]) / (2.0 * dx);

      // TODO : delete below once confirmed. Old way to get gradient
      //        that might have indexing issue
      //  Ratio of consecutive gradients (upwind-side / downwind-side)
      // const double du_up = (u[2] - u[1]);
      // const double du_down = (u[3] - u[2]);
      // const double r = du_down / (du_up + ((du_up >= 0.0) ? eps : -eps));

      // Slopes used by the left-biased stencil:
      //
      //   u[i-2] ---- u[i-1] ---- u[i]
      //     u[0]        u[1]      u[2]
      //
      // The second-order correction is proportional to
      //   (u[i] - u[i-1]) - (u[i-1] - u[i-2]).
      //
      // Therefore the limiter ratio should compare the previous upwind
      // slope to the current upwind slope.
      const double du_curr = u[2] - u[1]; // u_i     - u_{i-1}
      const double du_prev = u[1] - u[0]; // u_{i-1} - u_{i-2}
      const double r = du_prev / (du_curr + ((du_curr >= 0.0) ? eps : -eps));

      const double limiter = limiter_minmod(r);
      return d1 + limiter * (d2 - d1);
    }
    else {
      // a < 0 : upwind is to the right
      const double d1 = (u[3] - u[2]) / dx;

      // 2nd order upwind (biased right) derivative at i
      const double d2 = (-3.0 * u[2] + 4.0 * u[3] - u[4]) / (2.0 * dx);

      // const double du_up = (u[3] - u[2]);
      // const double du_down = (u[2] - u[1]);
      // const double r = du_down / (du_up + ((du_up >= 0.0) ? eps : -eps));

      // Slopes used by the right-biased stencil:
      //
      //   u[i] ---- u[i+1] ---- u[i+2]
      //   u[2]       u[3]        u[4]
      //
      // The second-order correction is proportional to
      //   (u[i+1] - u[i]) - (u[i+2] - u[i+1]).
      //
      // Therefore compare the next upwind slope to the current upwind slope.
      const double du_curr = u[3] - u[2]; // u_{i+1} - u_i
      const double du_next = u[4] - u[3]; // u_{i+2} - u_{i+1}
      const double r = du_next / (du_curr + ((du_curr >= 0.0) ? eps : -eps));

      const double limiter = limiter_minmod(r);
      return d1 + limiter * (d2 - d1);
    }
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

  // Compute deltas for each dimension
  std::array<double, D> div_d;
  if constexpr(D >= 1)
    div_d[0] = 1.0 / (2.0 * m.template delta<ax::x>());
  if constexpr(D >= 2)
    div_d[1] = 1.0 / (2.0 * m.template delta<ax::y>());
  if constexpr(D >= 3)
    div_d[2] = 1.0 / (2.0 * m.template delta<ax::z>());

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dp(i).x() = (psi(i + 1) - psi(i - 1)) * div_d[0];
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      dp(i, j).x() = (psi(i + 1, j) - psi(i - 1, j)) * div_d[0];
      dp(i, j).y() = (psi(i, j + 1) - psi(i, j - 1)) * div_d[1];
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      dp(i, j, k).x() = (psi(i + 1, j, k) - psi(i - 1, j, k)) * div_d[0];
      dp(i, j, k).y() = (psi(i, j + 1, k) - psi(i, j - 1, k)) * div_d[1];
      dp(i, j, k).z() = (psi(i, j, k + 1) - psi(i, j, k - 1)) * div_d[2];
    };
  }
}

struct unit_normals_task : annotation_task_region {
  inline static const std::string name{"unit_normals"};
};

// Unit normal computation: n = \nabla psi / |\nabla psi|
template<std::size_t D>
void
unit_normals(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<vec<D>>::template accessor<ro, ro> dp_a,
  typename field<vec<D>>::template accessor<wo, na> normal_a,
  typename faces_vec<D>::template accessor<wo, na> normal_f_a,
  std::size_t face_axis) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto dp = m.template mdcolex<is::cells>(dp_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);
  auto [n_right, n_left] = faces_vec<D>::mdcolex(m, normal_f_a);

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

      // Right face of cell i: x_{i+1/2}
      n_right(i) = normalize(0.5 * (dp(i) + dp(i + 1)));
      // Left face of cell i: x_{i-1/2}
      n_left(i) = normalize(0.5 * (dp(i - 1) + dp(i)));
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(dp,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      normal(i, j) = dp(i, j) / (dp(i, j).norm() + eps);

      if(face_axis == 0) {
        // x faces: i +/- 1/2
        n_right(i, j) = normalize(0.5 * (dp(i, j) + dp(i + 1, j)));
        n_left(i, j) = normalize(0.5 * (dp(i - 1, j) + dp(i, j)));
      }
      else {
        // y faces: j +/- 1/2
        n_right(i, j) = normalize(0.5 * (dp(i, j) + dp(i, j + 1)));
        n_left(i, j) = normalize(0.5 * (dp(i, j - 1) + dp(i, j)));
      }
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

      if(face_axis == 0) {
        n_right(i, j, k) = normalize(0.5 * (dp(i, j, k) + dp(i + 1, j, k)));
        n_left(i, j, k) = normalize(0.5 * (dp(i - 1, j, k) + dp(i, j, k)));
      }
      else if(face_axis == 1) {
        n_right(i, j, k) = normalize(0.5 * (dp(i, j, k) + dp(i, j + 1, k)));
        n_left(i, j, k) = normalize(0.5 * (dp(i, j - 1, k) + dp(i, j, k)));
      }
      else {
        n_right(i, j, k) = normalize(0.5 * (dp(i, j, k) + dp(i, j, k + 1)));
        n_left(i, j, k) = normalize(0.5 * (dp(i, j, k - 1) + dp(i, j, k)));
      }
    };
  }
}

// Compute directional derivative q_n = n · \nabla q, in the known region (psi
// <= 0) This corresponds to Eq. (4) in Aslam's paper.
template<std::size_t D>
void
compute_directional_derivative(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> q_a,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    qn_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto q = m.template mdcolex<is::cells>(q_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);
  auto qn = m.template mdcolex<is::cells>(qn_a);

  const double dx = m.template delta<ax::x>();

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(psi(i) > 0.0) {
        qn(i) = 0.0; // Only compute in known region (psi <= 0)
      }
      else {
        for(std::size_t qid = 0; qid < q(i).size(); ++qid) {
          // q_n = n · \nabla q,  using central differences
          qn(i)[qid] =
            normal(i).x() * (q(i + 1)[qid] - q(i - 1)[qid]) / (2.0 * dx);
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
      if(psi(i, j) > 0.0) {
        qn(i, j) = 0.0;
      }
      else {
        for(std::size_t qid = 0; qid < q(i, j).size(); ++qid) {
          // q_n = n · \nabla q = nx * \del q / \del x + ny * \del q / \del y
          qn(i, j)[qid] = normal(i, j).x() *
                            (q(i + 1, j)[qid] - q(i - 1, j)[qid]) / (2.0 * dx) +
                          normal(i, j).y() *
                            (q(i, j + 1)[qid] - q(i, j - 1)[qid]) / (2.0 * dy);
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
      if(psi(i, j, k) > 0.0) {
        qn(i, j, k) = 0.0;
      }
      else {
        for(std::size_t qid = 0; qid < q(i, j, k).size(); ++qid) {
          // q_n = n · \nabla q = nx * \del q / \del x + ny * \del q / \del y +
          // nz * \del q / \del z
          qn(i, j, k)[qid] =
            normal(i, j, k).x() * (q(i + 1, j, k)[qid] - q(i - 1, j, k)[qid]) /
              (2.0 * dx) +
            normal(i, j, k).y() * (q(i, j + 1, k)[qid] - q(i, j - 1, k)[qid]) /
              (2.0 * dy) +
            normal(i, j, k).z() * (q(i, j, k + 1)[qid] - q(i, j, k - 1)[qid]) /
              (2.0 * dz);
        }
      }
    };
  }
}

// Compute second directional derivative q_nn = n̂ · ∇(n̂ · ∇q) in the known
// region (ψ ≤ 0) This corresponds to Eq. (7) in Aslam's paper: u_nn = n̂ · ∇(n̂ ·
// ∇u)
template<std::size_t D>
void
compute_second_directional_derivative(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> q_a,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    qnn_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto q = m.template mdcolex<is::cells>(q_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);
  auto qnn = m.template mdcolex<is::cells>(qnn_a);

  const double dx = m.template delta<ax::x>();

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(psi(i) > 0.0) {
        qnn(i) = 0.0; // Only compute in known region (ψ ≤ 0)
      }
      else {
        for(std::size_t qid = 0; qid < q(i).size(); ++qid) {
          // First compute q_n = n̂ · ∇q at i-1, i, i+1
          double qn_im1 =
            normal(i - 1).x() * (q(i)[qid] - q(i - 2)[qid]) / (2.0 * dx);
          double qn_i =
            normal(i).x() * (q(i + 1)[qid] - q(i - 1)[qid]) / (2.0 * dx);
          double qn_ip1 =
            normal(i + 1).x() * (q(i + 2)[qid] - q(i)[qid]) / (2.0 * dx);

          // Then compute q_nn = n̂ · ∇(q_n) using central differences
          qnn(i)[qid] = normal(i).x() * (qn_ip1 - qn_im1) / (2.0 * dx);
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
      if(psi(i, j) > 0.0) {
        qnn(i, j) = 0.0;
      }
      else {
        for(std::size_t qid = 0; qid < q(i, j).size(); ++qid) {
          // Compute q_n at neighboring cells
          // For x-direction neighbors
          double qn_x_im1 = normal(i - 1, j).x() *
                              (q(i, j)[qid] - q(i - 2, j)[qid]) / (2.0 * dx) +
                            normal(i - 1, j).y() *
                              (q(i - 1, j + 1)[qid] - q(i - 1, j - 1)[qid]) /
                              (2.0 * dy);
          double qn_x_i = normal(i, j).x() *
                            (q(i + 1, j)[qid] - q(i - 1, j)[qid]) / (2.0 * dx) +
                          normal(i, j).y() *
                            (q(i, j + 1)[qid] - q(i, j - 1)[qid]) / (2.0 * dy);
          double qn_x_ip1 = normal(i + 1, j).x() *
                              (q(i + 2, j)[qid] - q(i, j)[qid]) / (2.0 * dx) +
                            normal(i + 1, j).y() *
                              (q(i + 1, j + 1)[qid] - q(i + 1, j - 1)[qid]) /
                              (2.0 * dy);

          // For y-direction neighbors
          double qn_y_jm1 = normal(i, j - 1).x() *
                              (q(i + 1, j - 1)[qid] - q(i - 1, j - 1)[qid]) /
                              (2.0 * dx) +
                            normal(i, j - 1).y() *
                              (q(i, j)[qid] - q(i, j - 2)[qid]) / (2.0 * dy);
          double qn_y_j = normal(i, j).x() *
                            (q(i + 1, j)[qid] - q(i - 1, j)[qid]) / (2.0 * dx) +
                          normal(i, j).y() *
                            (q(i, j + 1)[qid] - q(i, j - 1)[qid]) / (2.0 * dy);
          double qn_y_jp1 = normal(i, j + 1).x() *
                              (q(i + 1, j + 1)[qid] - q(i - 1, j + 1)[qid]) /
                              (2.0 * dx) +
                            normal(i, j + 1).y() *
                              (q(i, j + 2)[qid] - q(i, j)[qid]) / (2.0 * dy);

          // q_nn = n̂ · ∇(q_n) = nx * ∂q_n/∂x + ny * ∂q_n/∂y
          qnn(i, j)[qid] =
            normal(i, j).x() * (qn_x_ip1 - qn_x_im1) / (2.0 * dx) +
            normal(i, j).y() * (qn_y_jp1 - qn_y_jm1) / (2.0 * dy);
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
      if(psi(i, j, k) > 0.0) {
        qnn(i, j, k) = 0.0;
      }
      else {
        for(std::size_t qid = 0; qid < q(i, j, k).size(); ++qid) {
          // Compute q_n at neighboring cells in x-direction
          double qn_x_im1 =
            normal(i - 1, j, k).x() * (q(i, j, k)[qid] - q(i - 2, j, k)[qid]) /
              (2.0 * dx) +
            normal(i - 1, j, k).y() *
              (q(i - 1, j + 1, k)[qid] - q(i - 1, j - 1, k)[qid]) / (2.0 * dy) +
            normal(i - 1, j, k).z() *
              (q(i - 1, j, k + 1)[qid] - q(i - 1, j, k - 1)[qid]) / (2.0 * dz);
          double qn_x_i =
            normal(i, j, k).x() * (q(i + 1, j, k)[qid] - q(i - 1, j, k)[qid]) /
              (2.0 * dx) +
            normal(i, j, k).y() * (q(i, j + 1, k)[qid] - q(i, j - 1, k)[qid]) /
              (2.0 * dy) +
            normal(i, j, k).z() * (q(i, j, k + 1)[qid] - q(i, j, k - 1)[qid]) /
              (2.0 * dz);
          double qn_x_ip1 =
            normal(i + 1, j, k).x() * (q(i + 2, j, k)[qid] - q(i, j, k)[qid]) /
              (2.0 * dx) +
            normal(i + 1, j, k).y() *
              (q(i + 1, j + 1, k)[qid] - q(i + 1, j - 1, k)[qid]) / (2.0 * dy) +
            normal(i + 1, j, k).z() *
              (q(i + 1, j, k + 1)[qid] - q(i + 1, j, k - 1)[qid]) / (2.0 * dz);

          // Compute q_n at neighboring cells in y-direction
          double qn_y_jm1 =
            normal(i, j - 1, k).x() *
              (q(i + 1, j - 1, k)[qid] - q(i - 1, j - 1, k)[qid]) / (2.0 * dx) +
            normal(i, j - 1, k).y() * (q(i, j, k)[qid] - q(i, j - 2, k)[qid]) /
              (2.0 * dy) +
            normal(i, j - 1, k).z() *
              (q(i, j - 1, k + 1)[qid] - q(i, j - 1, k - 1)[qid]) / (2.0 * dz);
          double qn_y_j =
            normal(i, j, k).x() * (q(i + 1, j, k)[qid] - q(i - 1, j, k)[qid]) /
              (2.0 * dx) +
            normal(i, j, k).y() * (q(i, j + 1, k)[qid] - q(i, j - 1, k)[qid]) /
              (2.0 * dy) +
            normal(i, j, k).z() * (q(i, j, k + 1)[qid] - q(i, j, k - 1)[qid]) /
              (2.0 * dz);
          double qn_y_jp1 =
            normal(i, j + 1, k).x() *
              (q(i + 1, j + 1, k)[qid] - q(i - 1, j + 1, k)[qid]) / (2.0 * dx) +
            normal(i, j + 1, k).y() * (q(i, j + 2, k)[qid] - q(i, j, k)[qid]) /
              (2.0 * dy) +
            normal(i, j + 1, k).z() *
              (q(i, j + 1, k + 1)[qid] - q(i, j + 1, k - 1)[qid]) / (2.0 * dz);

          // Compute q_n at neighboring cells in z-direction
          double qn_z_km1 =
            normal(i, j, k - 1).x() *
              (q(i + 1, j, k - 1)[qid] - q(i - 1, j, k - 1)[qid]) / (2.0 * dx) +
            normal(i, j, k - 1).y() *
              (q(i, j + 1, k - 1)[qid] - q(i, j - 1, k - 1)[qid]) / (2.0 * dy) +
            normal(i, j, k - 1).z() * (q(i, j, k)[qid] - q(i, j, k - 2)[qid]) /
              (2.0 * dz);
          double qn_z_k =
            normal(i, j, k).x() * (q(i + 1, j, k)[qid] - q(i - 1, j, k)[qid]) /
              (2.0 * dx) +
            normal(i, j, k).y() * (q(i, j + 1, k)[qid] - q(i, j - 1, k)[qid]) /
              (2.0 * dy) +
            normal(i, j, k).z() * (q(i, j, k + 1)[qid] - q(i, j, k - 1)[qid]) /
              (2.0 * dz);
          double qn_z_kp1 =
            normal(i, j, k + 1).x() *
              (q(i + 1, j, k + 1)[qid] - q(i - 1, j, k + 1)[qid]) / (2.0 * dx) +
            normal(i, j, k + 1).y() *
              (q(i, j + 1, k + 1)[qid] - q(i, j - 1, k + 1)[qid]) / (2.0 * dy) +
            normal(i, j, k + 1).z() * (q(i, j, k + 2)[qid] - q(i, j, k)[qid]) /
              (2.0 * dz);

          // q_nn = n̂ · ∇(q_n) = nx * ∂q_n/∂x + ny * ∂q_n/∂y + nz * ∂q_n/∂z
          qnn(i, j, k)[qid] =
            normal(i, j, k).x() * (qn_x_ip1 - qn_x_im1) / (2.0 * dx) +
            normal(i, j, k).y() * (qn_y_jp1 - qn_y_jm1) / (2.0 * dy) +
            normal(i, j, k).z() * (qn_z_kp1 - qn_z_km1) / (2.0 * dz);
        }
      }
    };
  }
}

// -------------------------------------------------------------------
// RHS: R(u) = - H(psi) * (n · grad(u))  in psi>0, else 0.
// grad(u) computed with 2nd-order upwind + MinMod for each component.
// -------------------------------------------------------------------

struct rhs_constant_extrapolation_minmod_task : annotation_task_region {
  inline static const std::string name{"rhs_constant_extrapolation_minmod"};
};

template<std::size_t D>
void
rhs_constant_extrapolation_minmod(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro>
    q_scalar_a,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> normal_a,
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
        const double ax = normal(i).x();
        for(std::size_t qid = 0; qid < q(i).size(); ++qid) {
          // clang-format off
          rhs(i)[qid] = -ax * gradu_1d<5>({q(i - 2)[qid], q(i - 1)[qid],
              q(i)[qid], q(i + 1)[qid], q(i + 2)[qid]}, ax, dx);
          // clang-format on
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
        const double ax = normal(i, j).x();
        const double ay = normal(i, j).y();

        for(std::size_t qid = 0; qid < q(i, j).size(); ++qid) {
          // clang-format off
          rhs(i, j)[qid] = -(
            ax * gradu_1d<5>({q(i - 2, j)[qid], q(i - 1, j)[qid], q(i, j)[qid],
                               q(i + 1, j)[qid], q(i + 2, j)[qid]}, ax, dx) +
            ay * gradu_1d<5>({q(i, j - 2)[qid], q(i, j - 1)[qid], q(i, j)[qid],
                               q(i, j + 1)[qid], q(i, j + 2)[qid]}, ay, dy));
          // clang-format on
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
        const double ax = normal(i, j, k).x();
        const double ay = normal(i, j, k).y();
        const double az = normal(i, j, k).z();

        for(std::size_t qid = 0; qid < q(i, j, k).size(); ++qid) {
          // clang-format off
          rhs(i, j, k)[qid] = -(
            ax * gradu_1d<5>({q(i - 2, j, k)[qid], q(i - 1, j, k)[qid],
                               q(i, j, k)[qid], q(i + 1, j, k)[qid],
                               q(i + 2, j, k)[qid]}, ax, dx) +
            ay * gradu_1d<5>({q(i, j - 2, k)[qid], q(i, j - 1, k)[qid],
                               q(i, j, k)[qid], q(i, j + 1, k)[qid],
                               q(i, j + 2, k)[qid]}, ay, dy) +
            az * gradu_1d<5>({q(i, j, k - 2)[qid], q(i, j, k - 1)[qid],
                               q(i, j, k)[qid], q(i, j, k + 1)[qid],
                               q(i, j, k + 2)[qid]}, az, dz));
          // clang-format on
        }
      }
    };
  }
}

// -------------------------------------------------------------------
// RHS: R(u_n) = - H(psi + \delta x) * (n · grad(u_n))  in psi>0, else 0.
// grad(u_n) computed with 2nd-order upwind + MinMod for each component.
// -------------------------------------------------------------------

// RHS for linear extrapolation: ∂u/∂t = -H(ψ)(n̂ · ∇u - u_n)
// Corresponds to Eq. (6) in Aslam's paper
template<std::size_t D>
void
rhs_linear_extrapolation_minmod(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> q_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> qn_a,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    rhs_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto q = m.template mdcolex<is::cells>(q_a);
  auto qn = m.template mdcolex<is::cells>(qn_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto normal = m.template mdcolex<is::cells>(normal_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);

  const double dx = m.template delta<ax::x>();

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(psi(i) <= 0.0) {
        rhs(i) = 0.0; // H(ψ) = 0 in known region
      }
      else {
        const double ax = normal(i).x();
        for(std::size_t qid = 0; qid < q(i).size(); ++qid) {
          double grad_u = gradu_1d<5>({q(i - 2)[qid],
                                        q(i - 1)[qid],
                                        q(i)[qid],
                                        q(i + 1)[qid],
                                        q(i + 2)[qid]},
            ax,
            dx);

          // RHS = -H(ψ)(n̂ · ∇u - u_n) = -ax * grad_u + qn
          rhs(i)[qid] = -ax * grad_u + qn(i)[qid];
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
        const double ax = normal(i, j).x();
        const double ay = normal(i, j).y();

        for(std::size_t qid = 0; qid < q(i, j).size(); ++qid) {
          double grad_u_x = gradu_1d<5>({q(i - 2, j)[qid],
                                          q(i - 1, j)[qid],
                                          q(i, j)[qid],
                                          q(i + 1, j)[qid],
                                          q(i + 2, j)[qid]},
            ax,
            dx);
          double grad_u_y = gradu_1d<5>({q(i, j - 2)[qid],
                                          q(i, j - 1)[qid],
                                          q(i, j)[qid],
                                          q(i, j + 1)[qid],
                                          q(i, j + 2)[qid]},
            ay,
            dy);

          // RHS = -(ax * ∂q/∂x + ay * ∂q/∂y) + q_n
          rhs(i, j)[qid] = -(ax * grad_u_x + ay * grad_u_y) + qn(i, j)[qid];
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
        const double ax = normal(i, j, k).x();
        const double ay = normal(i, j, k).y();
        const double az = normal(i, j, k).z();

        for(std::size_t qid = 0; qid < q(i, j, k).size(); ++qid) {
          double grad_u_x = gradu_1d<5>({q(i - 2, j, k)[qid],
                                          q(i - 1, j, k)[qid],
                                          q(i, j, k)[qid],
                                          q(i + 1, j, k)[qid],
                                          q(i + 2, j, k)[qid]},
            ax,
            dx);
          double grad_u_y = gradu_1d<5>({q(i, j - 2, k)[qid],
                                          q(i, j - 1, k)[qid],
                                          q(i, j, k)[qid],
                                          q(i, j + 1, k)[qid],
                                          q(i, j + 2, k)[qid]},
            ay,
            dy);
          double grad_u_z = gradu_1d<5>({q(i, j, k - 2)[qid],
                                          q(i, j, k - 1)[qid],
                                          q(i, j, k)[qid],
                                          q(i, j, k + 1)[qid],
                                          q(i, j, k + 2)[qid]},
            az,
            dz);

          // RHS = -(ax * ∂q/∂x + ay * ∂q/∂y + az * ∂q/∂z) + q_n
          rhs(i, j, k)[qid] =
            -(ax * grad_u_x + ay * grad_u_y + az * grad_u_z) + qn(i, j, k)[qid];
        }
      }
    };
  }
}

// RHS for quadratic extrapolation of u_n: ∂u_n/∂t = -H(ψ)(n̂ · ∇u_n - u_nn)
// Corresponds to Eq. (9) in Aslam's paper
template<std::size_t D>
void
rhs_quadratic_extrapolation_qn_minmod(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> qn_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> qnn_a,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> normal_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    rhs_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto qn = m.template mdcolex<is::cells>(qn_a);
  auto qnn = m.template mdcolex<is::cells>(qnn_a);
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
        const double ax = normal(i).x();
        for(std::size_t qid = 0; qid < qn(i).size(); ++qid) {
          double grad_qn = gradu_1d<5>({qn(i - 2)[qid],
                                         qn(i - 1)[qid],
                                         qn(i)[qid],
                                         qn(i + 1)[qid],
                                         qn(i + 2)[qid]},
            ax,
            dx);

          // RHS = -H(ψ)(n̂ · ∇u_n - u_nn) = -ax * grad_qn + qnn
          rhs(i)[qid] = -ax * grad_qn + qnn(i)[qid];
        }
      }
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(qn,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    const double dy = m.template delta<ax::y>();

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      if(psi(i, j) <= 0.0) {
        rhs(i, j) = 0.0;
      }
      else {
        const double ax = normal(i, j).x();
        const double ay = normal(i, j).y();

        for(std::size_t qid = 0; qid < qn(i, j).size(); ++qid) {
          double grad_qn_x = gradu_1d<5>({qn(i - 2, j)[qid],
                                           qn(i - 1, j)[qid],
                                           qn(i, j)[qid],
                                           qn(i + 1, j)[qid],
                                           qn(i + 2, j)[qid]},
            ax,
            dx);
          double grad_qn_y = gradu_1d<5>({qn(i, j - 2)[qid],
                                           qn(i, j - 1)[qid],
                                           qn(i, j)[qid],
                                           qn(i, j + 1)[qid],
                                           qn(i, j + 2)[qid]},
            ay,
            dy);

          // RHS = -(ax * ∂q_n/∂x + ay * ∂q_n/∂y) + q_nn
          rhs(i, j)[qid] = -(ax * grad_qn_x + ay * grad_qn_y) + qnn(i, j)[qid];
        }
      }
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(qn,
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
        const double ax = normal(i, j, k).x();
        const double ay = normal(i, j, k).y();
        const double az = normal(i, j, k).z();

        for(std::size_t qid = 0; qid < qn(i, j, k).size(); ++qid) {
          double grad_qn_x = gradu_1d<5>({qn(i - 2, j, k)[qid],
                                           qn(i - 1, j, k)[qid],
                                           qn(i, j, k)[qid],
                                           qn(i + 1, j, k)[qid],
                                           qn(i + 2, j, k)[qid]},
            ax,
            dx);
          double grad_qn_y = gradu_1d<5>({qn(i, j - 2, k)[qid],
                                           qn(i, j - 1, k)[qid],
                                           qn(i, j, k)[qid],
                                           qn(i, j + 1, k)[qid],
                                           qn(i, j + 2, k)[qid]},
            ay,
            dy);
          double grad_qn_z = gradu_1d<5>({qn(i, j, k - 2)[qid],
                                           qn(i, j, k - 1)[qid],
                                           qn(i, j, k)[qid],
                                           qn(i, j, k + 1)[qid],
                                           qn(i, j, k + 2)[qid]},
            az,
            dz);

          // RHS = -(ax * ∂q_n/∂x + ay * ∂q_n/∂y + az * ∂q_n/∂z) + q_nn
          rhs(i, j, k)[qid] =
            -(ax * grad_qn_x + ay * grad_qn_y + az * grad_qn_z) +
            qnn(i, j, k)[qid];
        }
      }
    };
  }
}

// Pseudo-time step helper
template<std::size_t D>
inline double
min_mesh_spacing(typename mesh<D>::template accessor<ro> m) noexcept {
  if constexpr(D == 1) {
    return m.template delta<ax::x>();
  }
  else if constexpr(D == 2) {
    return std::min(m.template delta<ax::x>(), m.template delta<ax::y>());
  }
  else {
    return std::min({m.template delta<ax::x>(),
      m.template delta<ax::y>(),
      m.template delta<ax::z>()});
  }
}

template<std::size_t D>
FLECSI_INLINE_TARGET hard::tasks::util::extp_q<D>
q_stage1(const hard::tasks::util::extp_q<D> & q0,
  const hard::tasks::util::extp_q<D> & rhs,
  double dtau) noexcept {

  hard::tasks::util::extp_q<D> out{};

  out.rho() = q0.rho() + dtau * rhs.rho();

  for(std::size_t d = 0; d < D; ++d) {
    out.u()[d] = q0.u()[d] + dtau * rhs.u()[d];
  }

  out.p() = q0.p() + dtau * rhs.p();
  out.e() = q0.e() + dtau * rhs.e();
  out.n_comps() = q0.n_comps();

  return out;
}

template<std::size_t D>
FLECSI_INLINE_TARGET hard::tasks::util::extp_q<D>
q_stage2(const hard::tasks::util::extp_q<D> & q0,
  const hard::tasks::util::extp_q<D> & q1,
  const hard::tasks::util::extp_q<D> & rhs1,
  double dtau) noexcept {

  hard::tasks::util::extp_q<D> out{};

  out.rho() = 0.5 * (q0.rho() + q1.rho() + dtau * rhs1.rho());

  for(std::size_t d = 0; d < D; ++d) {
    out.u()[d] = 0.5 * (q0.u()[d] + q1.u()[d] + dtau * rhs1.u()[d]);
  }

  out.p() = 0.5 * (q0.p() + q1.p() + dtau * rhs1.p());
  out.e() = 0.5 * (q0.e() + q1.e() + dtau * rhs1.e());
  out.n_comps() = q0.n_comps();

  return out;
}

// PREVIOUS RK2, CHANGING IT TO THE CFL DEPENDENT ONE RIGHT BELOW
// RK2(Heun) advance for constant extrapolation
// STEP=1: q_stage = q + dt * rhs
// STEP=2: q_out   = 0.5*( q + q_stage + dt*rhs_stage )
template<std::size_t D>
void
advance_q_rk2_s1(flecsi::exec::accelerator s,
  typename field<double>::template accessor<ro, na> psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_a, // base q^n
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    rhs_a, // rhs at base or stage
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_new_a, // q^(1) for STEP=1, q^(n+1) for STEP=2
  typename single<double>::template accessor<ro> dt_a) noexcept {

  s.executor().forall(i, flecsi::util::iota_view({}, q_new_a.span().size())) {
    const double dt = *dt_a;
    if(psi_a(i) <= 0.0) {
      q_new_a(i) = q_known_a(i);
    }
    else {
      for(std::size_t qid = 0; qid < q_a(i).size(); ++qid) {
        q_new_a(i)[qid] = q_a(i)[qid] + dt * rhs_a(i)[qid];
      }
    }
  };
}

template<std::size_t D>
void
advance_q_rk2_s2(flecsi::exec::accelerator s,
  typename field<double>::template accessor<ro, na> psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_a, // base q^n
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    rhs_a, // rhs at base or stage
  typename field<hard::tasks::util::extp_q<D>>::template accessor<rw, na>
    q_stage_a, // q^(1) (only used in STEP=2)
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_new_a, // q^(1) for STEP=1, q^(n+1) for STEP=2
  typename single<double>::template accessor<ro> dt_a) noexcept {

  s.executor().forall(i, flecsi::util::iota_view({}, q_new_a.span().size())) {
    const double dt = *dt_a;
    if(psi_a(i) <= 0.0) {
      q_new_a(i) = q_known_a(i);
    }
    else {
      for(std::size_t qid = 0; qid < q_a(i).size(); ++qid) {
        q_new_a(i)[qid] =
          0.5 * (q_a(i)[qid] + q_stage_a(i)[qid] + dt * rhs_a(i)[qid]);
      }
    }
  };
}

struct advance_q_rk2_s1_cfl_task : annotation_task_region {
  inline static const std::string name{"advance_q_rk2_s1_cfl"};
};

template<std::size_t D>
void
advance_q_rk2_s1_cfl(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> extp_psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> q0_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> rhs_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na> q1_a,
  double extp_cfl) noexcept {

  const double dtau = extp_cfl * min_mesh_spacing<D>(m);

  s.executor().forall(i, flecsi::util::iota_view({}, q_known_a.span().size())) {
    q1_a(i) = (extp_psi_a(i) <= 0.0) ? q_known_a(i)
                                     : q_stage1<D>(q0_a(i), rhs_a(i), dtau);
  };
}

struct advance_q_rk2_s2_cfl_task : annotation_task_region {
  inline static const std::string name{"advance_q_rk2_s2_cfl"};
};

template<std::size_t D>
void
advance_q_rk2_s2_cfl(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> extp_psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> q0_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    rhs1_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> q1_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    qnew_a,
  double extp_cfl) noexcept {

  const double dtau = extp_cfl * min_mesh_spacing<D>(m);

  s.executor().forall(i, flecsi::util::iota_view({}, q_known_a.span().size())) {
    qnew_a(i) = (extp_psi_a(i) <= 0.0)
                  ? q_known_a(i)
                  : q_stage2<D>(q0_a(i), q1_a(i), rhs1_a(i), dtau);
  };
}

struct copy_extp_q_task : annotation_task_region {
  inline static const std::string name{"copy_extp_q"};
};

template<std::size_t D>
void
copy_extp_q(flecsi::exec::accelerator s,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro> src_a,
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
//   signed extrapolation level set:
//     extp_psi <= 0 in source_material
//     extp_psi >  0 outside source_material
//
// q_scalar/q_scalar_known:
//   packed current primitive state. q_scalar is the extrapolation work field;
//   q_scalar_known is used by RK to keep the known region fixed.
template<std::size_t D>
void
prepare_extp_for_material(flecsi::exec::accelerator s,
  typename field<double>::template accessor<ro, na> ls_psi_a,
  typename field<int>::template accessor<ro, na> material_a,
  typename field<double>::template accessor<ro, na> mass_density_a,
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  typename field<double>::template accessor<ro, na> pressure_a,
  typename field<double>::template accessor<ro, na> sie_a,
  typename field<double>::template accessor<wo, na> extp_psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na> q_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_known_a,
  int source_material) noexcept {

  s.executor().forall(
    i, flecsi::util::iota_view({}, extp_psi_a.span().size())) {
    extp_psi_a(i) =
      (material_a(i) == source_material) ? ls_psi_a(i) : -ls_psi_a(i);

    q_a(i).rho() = mass_density_a(i);
    q_a(i).u() = velocity_a(i);
    q_a(i).p() = pressure_a(i);
    q_a(i).e() = sie_a(i);
    q_a(i).n_comps() = static_cast<double>(hard::tasks::util::extp_q<D>::N_q);

    q_known_a(i) = q_a(i);
  };
}

} // namespace hard::tasks::extrapolation

#endif