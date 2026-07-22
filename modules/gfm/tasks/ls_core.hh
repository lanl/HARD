#ifndef HARD_MODULES_GFM_TASKS_LS_CORE_HH
#define HARD_MODULES_GFM_TASKS_LS_CORE_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/mesh.hh"
#include "../state.hh"

using namespace hard::tasks::util;

namespace hard::tasks::level_set {

struct freeze_material_task : annotation_task_region {
  inline static const std::string name{"freeze_material"};
};

struct rhs_upwind_task : annotation_task_region {
  inline static const std::string name{"rhs_upwind"};
};

// Freeze material labels at the start of each step
template<std::size_t D>
void
freeze_material(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, na> material_a,
  typename field<int>::template accessor<wo, na> frozen_material_a) noexcept {

  s.executor().forall(
    i, flecsi::util::iota_view({}, frozen_material_a.span().size())) {
    frozen_material_a(i) = material_a(i);
  };
}

namespace hj_weno5m_stencil {

// This scheme is used in the level set advection. It was taken from:
// "Mapped weighted essentially non-oscillatory schemes: Achieving optimal order
// near critical points" by Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
// Journal of Computational Physics, 2005.

// This scheme is used in the paper:
// "An efficient level set method for tracking many materials"
// by Tariq D. Aslam, Eduardo Lozano, Journal of Computational Physics, 2026.

// This method is needed for better accuracy in the level set advection.

template<typename T>
FLECSI_INLINE_TARGET constexpr T
sq(const T x) noexcept {
  return x * x;
}

template<typename T>
FLECSI_INLINE_TARGET T
weno5m_map_weight(const T omega, const T d) noexcept {
  // Henrick-Aslam-Powers mapped WENO5 weight.
  // Fixed points at omega = 0, d, 1.
  const T num = omega * (d + d * d - T(3) * d * omega + omega * omega);
  const T den = d * d + omega * (T(1) - T(2) * d);
  return num / den;
}

template<typename T>
FLECSI_INLINE_TARGET T
weno5m_left_from_differences(const std::array<T, 5> & v, const T eps) noexcept {
  // v = {v_{i-2}, v_{i-1}, v_i, v_{i+1}, v_{i+2}}
  // where v_j = (psi_{j+1} - psi_j) / dx

  constexpr T d0 = T(1) / T(10);
  constexpr T d1 = T(6) / T(10);
  constexpr T d2 = T(3) / T(10);

  // Candidate polynomials
  const T q0 =
    (T(1) / T(3)) * v[0] - (T(7) / T(6)) * v[1] + (T(11) / T(6)) * v[2];
  const T q1 =
    -(T(1) / T(6)) * v[1] + (T(5) / T(6)) * v[2] + (T(1) / T(3)) * v[3];
  const T q2 =
    (T(1) / T(3)) * v[2] + (T(5) / T(6)) * v[3] - (T(1) / T(6)) * v[4];

  // Smoothness indicators
  const T beta0 = (T(13) / T(12)) * sq(v[0] - T(2) * v[1] + v[2]) +
                  (T(1) / T(4)) * sq(v[0] - T(4) * v[1] + T(3) * v[2]);

  const T beta1 = (T(13) / T(12)) * sq(v[1] - T(2) * v[2] + v[3]) +
                  (T(1) / T(4)) * sq(v[1] - v[3]);

  const T beta2 = (T(13) / T(12)) * sq(v[2] - T(2) * v[3] + v[4]) +
                  (T(1) / T(4)) * sq(T(3) * v[2] - T(4) * v[3] + v[4]);

  // Jiang-Shu nonlinear weights
  const T a0 = d0 / sq(eps + beta0);
  const T a1 = d1 / sq(eps + beta1);
  const T a2 = d2 / sq(eps + beta2);

  const T asum = a0 + a1 + a2;

  const T w0 = a0 / asum;
  const T w1 = a1 / asum;
  const T w2 = a2 / asum;

  // WENO5-M mapped weights
  const T wm0 = weno5m_map_weight(w0, d0);
  const T wm1 = weno5m_map_weight(w1, d1);
  const T wm2 = weno5m_map_weight(w2, d2);

  const T wmsum = wm0 + wm1 + wm2;

  const T wt0 = wm0 / wmsum;
  const T wt1 = wm1 / wmsum;
  const T wt2 = wm2 / wmsum;

  return wt0 * q0 + wt1 * q1 + wt2 * q2;
}

template<typename T>
FLECSI_INLINE_TARGET T
weno5m_right_from_differences(const std::array<T, 5> & v,
  const T eps) noexcept {
  // Mirror of left-biased reconstruction.
  // If v corresponds to {v_{i-3}, v_{i-2}, v_{i-1}, v_i, v_{i+1}},
  // this returns the right-biased derivative at i.
  return -weno5m_left_from_differences<T>(
    {-v[4], -v[3], -v[2], -v[1], -v[0]}, eps);
}

template<typename T>
FLECSI_INLINE_TARGET std::pair<T, T>
weno5m_hj_derivatives_from_points(const std::array<T, 7> & p,
  const T dx) noexcept {
  // p = {psi_{i-3}, psi_{i-2}, psi_{i-1}, psi_i,
  //      psi_{i+1}, psi_{i+2}, psi_{i+3}}
  // Build consecutive first differences, then reconstruct one-sided
  // Hamilton-Jacobi derivatives at cell i.

  const std::array<T, 6> d = {
    (p[1] - p[0]) / dx,
    (p[2] - p[1]) / dx,
    (p[3] - p[2]) / dx,
    (p[4] - p[3]) / dx,
    (p[5] - p[4]) / dx,
    (p[6] - p[5]) / dx,
  };

  const T eps = T(1e-8);

  // Left-biased derivative phi_x^- at i uses differences spanning
  // {i-5/2, i-3/2, i-1/2, i+1/2, i+3/2} -> d[0]..d[4].
  const T dminus =
    weno5m_left_from_differences<T>({d[0], d[1], d[2], d[3], d[4]}, eps);

  // Right-biased derivative phi_x^+ at i uses the mirrored stencil
  // spanning {i-3/2, i-1/2, i+1/2, i+3/2, i+5/2},
  // which corresponds to d[1]..d[5].
  const T dplus =
    weno5m_right_from_differences<T>({d[1], d[2], d[3], d[4], d[5]}, eps);

  return {dminus, dplus};
}

} // namespace hj_weno5m_stencil

// upwind HJ RHS using material-signed psi* values so that neighbor
// differences are taken within the center cell's frozen material frame.
template<std::size_t D>
void
rhs_upwind(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  typename field<int>::template accessor<ro, ro> frozen_material_a,
  typename field<double>::template accessor<wo, na> rhs_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto vel = m.template mdcolex<is::cells>(velocity_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);

  if constexpr(D == 1) {
    const double dx = m.template delta<ax::x>();

    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const int center_material = fmat(i);
      const double ux = vel(i).x();

      const auto psi_star = [&](int ii) {
        const double sign = (fmat(ii) == center_material) ? 1.0 : -1.0;
        return psi(ii) * sign;
      };

      // 7-point stencil in the center cell's frozen-material frame:
      // {psi*_{i-3}, ..., psi*_{i+3}}
      const std::array<double, 7> p = {
        psi_star(i - 3),
        psi_star(i - 2),
        psi_star(i - 1),
        psi_star(i),
        psi_star(i + 1),
        psi_star(i + 2),
        psi_star(i + 3),
      };

      const auto [psi_x_minus, psi_x_plus] =
        hj_weno5m_stencil::weno5m_hj_derivatives_from_points(p, dx);

      const double psi_x_upwind = (ux > 0.0) ? psi_x_minus : psi_x_plus;
      rhs(i) = -(ux * psi_x_upwind);
    };
  }
  else if constexpr(D == 2) {
    const double dx = m.template delta<ax::x>();
    const double dy = m.template delta<ax::y>();

    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      const int center_material = fmat(i, j);
      const auto u = vel(i, j);

      const auto psi_star = [&](int ii, int jj) {
        const double sign = (fmat(ii, jj) == center_material) ? 1.0 : -1.0;
        return psi(ii, jj) * sign;
      };

      const std::array<double, 7> px = {
        psi_star(i - 3, j),
        psi_star(i - 2, j),
        psi_star(i - 1, j),
        psi_star(i, j),
        psi_star(i + 1, j),
        psi_star(i + 2, j),
        psi_star(i + 3, j),
      };

      const std::array<double, 7> py = {
        psi_star(i, j - 3),
        psi_star(i, j - 2),
        psi_star(i, j - 1),
        psi_star(i, j),
        psi_star(i, j + 1),
        psi_star(i, j + 2),
        psi_star(i, j + 3),
      };

      const auto [psi_x_minus, psi_x_plus] =
        hj_weno5m_stencil::weno5m_hj_derivatives_from_points(px, dx);
      const auto [psi_y_minus, psi_y_plus] =
        hj_weno5m_stencil::weno5m_hj_derivatives_from_points(py, dy);

      const double psi_x_upwind = (u.x() > 0.0) ? psi_x_minus : psi_x_plus;
      const double psi_y_upwind = (u.y() > 0.0) ? psi_y_minus : psi_y_plus;

      rhs(i, j) = -(u.x() * psi_x_upwind + u.y() * psi_y_upwind);
    };
  }
  else { // D == 3
    const double dx = m.template delta<ax::x>();
    const double dy = m.template delta<ax::y>();
    const double dz = m.template delta<ax::z>();

    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      const int center_material = fmat(i, j, k);
      const auto u = vel(i, j, k);

      const auto psi_star = [&](int ii, int jj, int kk) {
        const double sign = (fmat(ii, jj, kk) == center_material) ? 1.0 : -1.0;
        return psi(ii, jj, kk) * sign;
      };

      const std::array<double, 7> px = {
        psi_star(i - 3, j, k),
        psi_star(i - 2, j, k),
        psi_star(i - 1, j, k),
        psi_star(i, j, k),
        psi_star(i + 1, j, k),
        psi_star(i + 2, j, k),
        psi_star(i + 3, j, k),
      };

      const std::array<double, 7> py = {
        psi_star(i, j - 3, k),
        psi_star(i, j - 2, k),
        psi_star(i, j - 1, k),
        psi_star(i, j, k),
        psi_star(i, j + 1, k),
        psi_star(i, j + 2, k),
        psi_star(i, j + 3, k),
      };

      const std::array<double, 7> pz = {
        psi_star(i, j, k - 3),
        psi_star(i, j, k - 2),
        psi_star(i, j, k - 1),
        psi_star(i, j, k),
        psi_star(i, j, k + 1),
        psi_star(i, j, k + 2),
        psi_star(i, j, k + 3),
      };

      const auto [psi_x_minus, psi_x_plus] =
        hj_weno5m_stencil::weno5m_hj_derivatives_from_points(px, dx);
      const auto [psi_y_minus, psi_y_plus] =
        hj_weno5m_stencil::weno5m_hj_derivatives_from_points(py, dy);
      const auto [psi_z_minus, psi_z_plus] =
        hj_weno5m_stencil::weno5m_hj_derivatives_from_points(pz, dz);

      const double psi_x_upwind = (u.x() > 0.0) ? psi_x_minus : psi_x_plus;
      const double psi_y_upwind = (u.y() > 0.0) ? psi_y_minus : psi_y_plus;
      const double psi_z_upwind = (u.z() > 0.0) ? psi_z_minus : psi_z_plus;

      rhs(i, j, k) =
        -(u.x() * psi_x_upwind + u.y() * psi_y_upwind + u.z() * psi_z_upwind);
    };
  }
}

struct advance_psi_s1_task : annotation_task_region {
  inline static const std::string name{"advance_psi_s1"};
};

// SSPRK3 stage update
template<std::size_t D, std::size_t STEP>
void
advance_psi_s1(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, na> psi_n0_a,
  typename field<double>::template accessor<ro, na> rhs_a,
  typename field<double>::template accessor<wo, na> psi_out_a,
  typename single<double>::template accessor<ro> dt_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi_n0 = m.template mdcolex<is::cells>(psi_n0_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);
  auto psi_out = m.template mdcolex<is::cells>(psi_out_a);

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      psi_out(i) = psi_n0(i) + *dt_a * rhs(i);
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      psi_out(i, j) = psi_n0(i, j) + *dt_a * rhs(i, j);
    };
  }
  else {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      psi_out(i, j, k) = psi_n0(i, j, k) + *dt_a * rhs(i, j, k);
    };
  }
}

struct advance_psi_task : annotation_task_region {
  inline static const std::string name{"advance_psi"};
};

template<std::size_t D, std::size_t STEP>
void
advance_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, na> psi_n0_a,
  typename field<double>::template accessor<ro, na> rhs_a,
  typename field<double>::template accessor<ro, na> psi_prev_a,
  typename field<double>::template accessor<wo, na> psi_out_a,
  typename single<double>::template accessor<ro> dt_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  static_assert(STEP >= 2 && STEP <= 3, "STEP must be 2 or 3 here");

  auto psi_n0 = m.template mdcolex<is::cells>(psi_n0_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);
  auto psi_prev = m.template mdcolex<is::cells>(psi_prev_a);
  auto psi_out = m.template mdcolex<is::cells>(psi_out_a);

  // const double dt = *dt_a;

  auto update_one = [=](double psi0, double psip, double k) {
    if constexpr(STEP == 2) {
      // psi^2 = 3/4 psi^n + 1/4 psi^1 + 1/4 dt L(psi^1, L^n)
      return 0.75 * psi0 + 0.25 * psip + 0.25 * (*dt_a) * k;
    }
    else {
      // psi^3 = 1/3 psi^n + 2/3 psi^2 + 2/3 dt L(psi^2, L^n)
      return (1.0 / 3.0) * psi0 + (2.0 / 3.0) * psip +
             (2.0 / 3.0) * (*dt_a) * k;
    }
  };

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      psi_out(i) = update_one(psi_n0(i), psi_prev(i), rhs(i));
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      psi_out(i, j) = update_one(psi_n0(i, j), psi_prev(i, j), rhs(i, j));
    };
  }
  else {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      psi_out(i, j, k) =
        update_one(psi_n0(i, j, k), psi_prev(i, j, k), rhs(i, j, k));
    };
  }
}

struct determine_psi_task : annotation_task_region {
  inline static const std::string name{"determine_psi"};
};

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
template<std::size_t D>
void
determine_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> frozen_psi_a, // psi3
  typename field<double>::template accessor<wo, na> psi_a, // psi^{n+1}
  typename field<int>::template accessor<ro, ro> frozen_material_a, // L^n
  typename field<int>::template accessor<wo, na> material_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto fpsi = m.template mdcolex<is::cells>(frozen_psi_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);
  auto mat = m.template mdcolex<is::cells>(material_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(fpsi(i) > 0.0) {
        double psimin = fpsi(i);
        int lmin = fmat(i);

        // Array of neighbor indices
        std::array<unsigned int, 2 * D> neighbors = {i - 1, i + 1};

        for(auto ni : neighbors) {
          if(fmat(ni) != fmat(i) && fpsi(ni) < psimin) {
            psimin = fpsi(ni);
            lmin = fmat(ni);
          }
        }

        psi(i) = (psimin >= 0.0 ? 1.0 : -1.0) * fpsi(i);
        mat(i) = lmin;
      }
      else {
        psi(i) = fpsi(i);
        mat(i) = fmat(i);
      }
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      if(fpsi(i, j) > 0.0) {
        double psimin = fpsi(i, j);
        int lmin = fmat(i, j);

        // Array of neighbor coordinate pairs
        std::array<std::tuple<unsigned int, unsigned int>, 2 * D> neighbors = {
          {{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}}};

        for(auto [ni, nj] : neighbors) {
          if(fmat(ni, nj) != fmat(i, j) && fpsi(ni, nj) < psimin) {
            psimin = fpsi(ni, nj);
            lmin = fmat(ni, nj);
          }
        }

        psi(i, j) = (psimin >= 0.0 ? 1.0 : -1.0) * fpsi(i, j);
        mat(i, j) = lmin;
      }
      else {
        psi(i, j) = fpsi(i, j);
        mat(i, j) = fmat(i, j);
      }
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      if(fpsi(i, j, k) > 0.0) {
        double psimin = fpsi(i, j, k);
        int lmin = fmat(i, j, k);

        // Array of neighbor coordinate triples
        std::array<std::tuple<unsigned int, unsigned int, unsigned int>, 2 * D>
          neighbors = {{{i - 1, j, k},
            {i + 1, j, k},
            {i, j - 1, k},
            {i, j + 1, k},
            {i, j, k - 1},
            {i, j, k + 1}}};

        for(auto [ni, nj, nk] : neighbors) {
          if(fmat(ni, nj, nk) != fmat(i, j, k) && fpsi(ni, nj, nk) < psimin) {
            psimin = fpsi(ni, nj, nk);
            lmin = fmat(ni, nj, nk);
          }
        }

        psi(i, j, k) = (psimin >= 0.0 ? 1.0 : -1.0) * fpsi(i, j, k);
        mat(i, j, k) = lmin;
      }
      else {
        psi(i, j, k) = fpsi(i, j, k);
        mat(i, j, k) = fmat(i, j, k);
      }
    };
  }
}

template<std::size_t D, typename eosType>
FLECSI_INLINE_TARGET bool
eos_pair_equivalent_at_cell(const eosType & eos_disp,
  int old_mat,
  int new_mat,
  double rho,
  const vec<D> & mom,
  double rhoE) noexcept {

  if(old_mat == new_mat) {
    return true;
  }

  if(old_mat < 0 || new_mat < 0 ||
     static_cast<std::size_t>(old_mat) >= eos_disp.size() ||
     static_cast<std::size_t>(new_mat) >= eos_disp.size()) {
    return false;
  }

  if(!std::isfinite(rho) || rho <= 0.0) {
    return false;
  }

  const vec<D> u = (1.0 / rho) * mom;
  const double e = (rhoE - 0.5 * rho * u.norm_squared()) / rho;

  if(!std::isfinite(e) || e <= 0.0) {
    return false;
  }

  const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
  const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

  const double p_old = eos_old.pRhoSie(rho, e);
  const double p_new = eos_new.pRhoSie(rho, e);
  const double c_old = eos_old.cRhoSie(rho, e);
  const double c_new = eos_new.cRhoSie(rho, e);

  constexpr double rtol = 1.0e-10;
  constexpr double atol = 1.0e-12;

  return std::abs(p_old - p_new) <=
           atol + rtol * std::max(std::abs(p_old), std::abs(p_new)) &&
         std::abs(c_old - c_new) <=
           atol + rtol * std::max(std::abs(c_old), std::abs(c_new));
}

struct promote_changed_material_cells_task : annotation_task_region {
  inline static const std::string name{"promote_changed_material_cells"};
};

template<std::size_t D>
void
promote_changed_material_cells(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, na> old_material_a,
  typename field<int>::template accessor<ro, na> new_material_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_ext_plus_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na>
    q_ext_minus_a,
  field<double>::accessor<rw, na> mass_density_a,
  typename field<vec<D>>::template accessor<rw, na> momentum_density_a,
  field<double>::accessor<rw, na> total_energy_density_a,
  typename field<vec<D>>::template accessor<rw, na> velocity_a,
  field<double>::accessor<rw, na> pressure_a,
  field<double>::accessor<rw, na> sie_a,
  field<double>::accessor<rw, na> sound_speed_a,
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {

  eos::eos_dispatch<flecsi::exec::accelerator> eos_disp(eos_mm);

  using common::tasks::utils::get_mdiota_policy;

  auto old_material = m.template mdcolex<is::cells>(old_material_a);
  auto new_material = m.template mdcolex<is::cells>(new_material_a);

  auto q_ext_plus = m.template mdcolex<is::cells>(q_ext_plus_a);
  auto q_ext_minus = m.template mdcolex<is::cells>(q_ext_minus_a);

  auto rho = m.template mdcolex<is::cells>(mass_density_a);
  auto mom = m.template mdcolex<is::cells>(momentum_density_a);
  auto rhoE = m.template mdcolex<is::cells>(total_energy_density_a);

  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
  auto sie = m.template mdcolex<is::cells>(sie_a);
  auto sound_speed = m.template mdcolex<is::cells>(sound_speed_a);

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      const int old_mat = old_material(i);
      const int new_mat = new_material(i);

      if(old_mat == new_mat) {
        return;
      }

      assert(
        new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_disp.size());

      const double r_old = rho(i);
      const vec<D> u_old = (1.0 / r_old) * mom(i);
      const double e_old =
        (rhoE(i) - 0.5 * r_old * u_old.norm_squared()) / r_old;

      const bool equivalent_eos = eos_pair_equivalent_at_cell<D>(
        eos_disp, old_mat, new_mat, r_old, mom(i), rhoE(i));

      if(equivalent_eos) {
        const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

        velocity(i) = u_old;
        sie(i) = e_old;
        pressure(i) = eos_new.pRhoSie(r_old, e_old);
        sound_speed(i) = eos_new.cRhoSie(r_old, e_old);

        return;
      }

      const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
      const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

      const double p_copy = eos_old.pRhoSie(r_old, e_old);
      const vec<D> u_copy = u_old;

      const auto & q = (new_mat == 0) ? q_ext_plus(i) : q_ext_minus(i);
      const double r_new = q.rho();

      assert(std::isfinite(r_new) && r_new > 0.0);
      assert(std::isfinite(p_copy) && p_copy > 0.0);

      const double e_new = hard::tasks::util::find_sie(eos_new, r_new, p_copy);
      const double c_new = eos_new.cRhoSie(r_new, e_new);

      rho(i) = r_new;
      velocity(i) = u_copy;
      pressure(i) = p_copy;
      sie(i) = e_new;
      sound_speed(i) = c_new;

      mom(i) = r_new * u_copy;
      rhoE(i) = r_new * e_new + 0.5 * r_new * u_copy.norm_squared();
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(rho,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      const int old_mat = old_material(i, j);
      const int new_mat = new_material(i, j);

      if(old_mat == new_mat) {
        return;
      }

      assert(
        new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_disp.size());

      const double r_old = rho(i, j);
      const vec<D> u_old = (1.0 / r_old) * mom(i, j);
      const double e_old =
        (rhoE(i, j) - 0.5 * r_old * u_old.norm_squared()) / r_old;

      const bool equivalent_eos = eos_pair_equivalent_at_cell<D>(
        eos_disp, old_mat, new_mat, r_old, mom(i, j), rhoE(i, j));

      if(equivalent_eos) {
        const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

        velocity(i, j) = u_old;
        sie(i, j) = e_old;
        pressure(i, j) = eos_new.pRhoSie(r_old, e_old);
        sound_speed(i, j) = eos_new.cRhoSie(r_old, e_old);

        return;
      }

      const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
      const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

      const double p_copy = eos_old.pRhoSie(r_old, e_old);
      const vec<D> u_copy = u_old;

      const auto & q = (new_mat == 0) ? q_ext_plus(i, j) : q_ext_minus(i, j);
      const double r_new = q.rho();

      assert(std::isfinite(r_new) && r_new > 0.0);
      assert(std::isfinite(p_copy) && p_copy > 0.0);

      const double e_new = hard::tasks::util::find_sie(eos_new, r_new, p_copy);
      const double c_new = eos_new.cRhoSie(r_new, e_new);

      rho(i, j) = r_new;
      velocity(i, j) = u_copy;
      pressure(i, j) = p_copy;
      sie(i, j) = e_new;
      sound_speed(i, j) = c_new;

      mom(i, j) = r_new * u_copy;
      rhoE(i, j) = r_new * e_new + 0.5 * r_new * u_copy.norm_squared();
    };
  }
  else {
    auto pol = get_mdiota_policy(rho,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      const int old_mat = old_material(i, j, k);
      const int new_mat = new_material(i, j, k);

      if(old_mat == new_mat) {
        return;
      }

      assert(
        new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_disp.size());

      const double r_old = rho(i, j, k);
      const vec<D> u_old = (1.0 / r_old) * mom(i, j, k);
      const double e_old =
        (rhoE(i, j, k) - 0.5 * r_old * u_old.norm_squared()) / r_old;

      const bool equivalent_eos = eos_pair_equivalent_at_cell<D>(
        eos_disp, old_mat, new_mat, r_old, mom(i, j, k), rhoE(i, j, k));

      if(equivalent_eos) {
        const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

        velocity(i, j, k) = u_old;
        sie(i, j, k) = e_old;
        pressure(i, j, k) = eos_new.pRhoSie(r_old, e_old);
        sound_speed(i, j, k) = eos_new.cRhoSie(r_old, e_old);

        return;
      }

      const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
      const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

      const double p_copy = eos_old.pRhoSie(r_old, e_old);
      const vec<D> u_copy = u_old;

      const auto & q =
        (new_mat == 0) ? q_ext_plus(i, j, k) : q_ext_minus(i, j, k);
      const double r_new = q.rho();

      assert(std::isfinite(r_new) && r_new > 0.0);
      assert(std::isfinite(p_copy) && p_copy > 0.0);

      const double e_new = hard::tasks::util::find_sie(eos_new, r_new, p_copy);
      const double c_new = eos_new.cRhoSie(r_new, e_new);

      rho(i, j, k) = r_new;
      velocity(i, j, k) = u_copy;
      pressure(i, j, k) = p_copy;
      sie(i, j, k) = e_new;
      sound_speed(i, j, k) = c_new;

      mom(i, j, k) = r_new * u_copy;
      rhoE(i, j, k) = r_new * e_new + 0.5 * r_new * u_copy.norm_squared();
    };
  }
}

// -----------------------------------------------------------------------------
// 1D diagnostic: print cells whose material changed during determine_psi.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
diagnose_material_changes_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, na> old_material_a,
  typename field<int>::template accessor<ro, na> new_material_a,
  field<double>::accessor<ro, na> mass_density_a,
  typename field<vec<D>>::template accessor<ro, na> momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  const std::vector<eos::eos_wrapper> & eos_mm,
  int max_print) noexcept {

  if constexpr(D == 1) {
    auto old_material = m.template mdcolex<is::cells>(old_material_a);
    auto new_material = m.template mdcolex<is::cells>(new_material_a);
    auto rho = m.template mdcolex<is::cells>(mass_density_a);
    auto mom = m.template mdcolex<is::cells>(momentum_density_a);
    auto rhoE = m.template mdcolex<is::cells>(total_energy_density_a);

    int nchanged = 0;
    int nprinted = 0;

    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const int old_mat = old_material(i);
      const int new_mat = new_material(i);

      if(old_mat == new_mat) {
        continue;
      }

      ++nchanged;

      const double r = rho(i);
      const double u = mom(i)[0] / r;
      const double e = (rhoE(i) - 0.5 * r * u * u) / r;

      double p_old = -1.0;
      double c_old = -1.0;
      double p_new = -1.0;
      double c_new = -1.0;

      if(old_mat >= 0 && static_cast<std::size_t>(old_mat) < eos_mm.size()) {
        const auto & eos_old = eos_mm[static_cast<std::size_t>(old_mat)];
        p_old = eos_old.pRhoSie(r, e);
        c_old = eos_old.cRhoSie(r, e);
      }

      if(new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_mm.size()) {
        const auto & eos_new = eos_mm[static_cast<std::size_t>(new_mat)];
        p_new = eos_new.pRhoSie(r, e);
        c_new = eos_new.cRhoSie(r, e);
      }

      if(nprinted < max_print) {
        std::cout << "[material changed] i=" << i << " old=" << old_mat
                  << " new=" << new_mat << " rho=" << r << " u=" << u
                  << " rhoE=" << rhoE(i) << " e_from_cons=" << e
                  << " p_oldEOS=" << p_old << " c_oldEOS=" << c_old
                  << " p_newEOS=" << p_new << " c_newEOS=" << c_new << "\n";

        ++nprinted;
      }
    }

    std::cout << "[material changed total] " << nchanged << "\n";
  }
}

struct diagnose_psi_range_1d_task : annotation_task_region {
  inline static const std::string name{"diagnose_psi_range_1d"};
};

template<std::size_t D>
void
diagnose_psi_range_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> psi_a,
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  const char * label) noexcept {

  if constexpr(D == 1) {
    auto psi = m.template mdcolex<is::cells>(psi_a);
    auto vel = m.template mdcolex<is::cells>(velocity_a);

    double min_psi = std::numeric_limits<double>::max();
    double max_psi = -std::numeric_limits<double>::max();
    double max_abs_psi = 0.0;
    double min_u = std::numeric_limits<double>::max();
    double max_u = -std::numeric_limits<double>::max();

    int bad = 0;

    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const double p = psi(i);
      const double u = vel(i)[0];

      min_psi = std::min(min_psi, p);
      max_psi = std::max(max_psi, p);
      max_abs_psi = std::max(max_abs_psi, std::abs(p));
      min_u = std::min(min_u, u);
      max_u = std::max(max_u, u);

      if(!std::isfinite(p) || std::abs(p) > 10.0) {
        ++bad;
        if(bad <= 10) {
          std::cout << "[bad psi] " << label << " i=" << i << " psi=" << p
                    << " u=" << u << "\n";
        }
      }
    }

    if(bad > 0) {
      std::cout << "[psi range] " << label << " min=" << min_psi
                << " max=" << max_psi << " maxabs=" << max_abs_psi
                << " u_min=" << min_u << " u_max=" << max_u
                << " bad_count=" << bad << "\n";
    }
  }
}

template<std::size_t D>
void
diagnose_material_cuts_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> material_a,
  field<double>::accessor<ro, ro> psi_a) noexcept {

  if constexpr(D == 1) {
    auto material = m.template mdcolex<is::cells>(material_a);
    auto psi = m.template mdcolex<is::cells>(psi_a);

    int ncuts = 0;

    for(auto i : m.template cells<ax::x, dm::corrector>()) {
      if(material(i - 1) != material(i)) {
        ++ncuts;

        const double xL = 0.5 * (m.template left<ax::x>(i - 1) +
                                  m.template right<ax::x>(i - 1));
        const double xR =
          0.5 * (m.template left<ax::x>(i) + m.template right<ax::x>(i));

        std::cout << "[material cut] i=" << i << " xL=" << xL << " xR=" << xR
                  << " mL=" << material(i - 1) << " mR=" << material(i)
                  << " psiL=" << psi(i - 1) << " psiR=" << psi(i) << "\n";
      }
    }

    std::cout << "[material cuts total] " << ncuts << "\n";
  }
}

} // namespace hard::tasks::level_set
#endif
