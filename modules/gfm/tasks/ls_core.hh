#ifndef HARD_MODULES_GFM_TASKS_LS_CORE_HH
#define HARD_MODULES_GFM_TASKS_LS_CORE_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/mesh.hh"
#include "../state.hh"

using namespace hard::tasks::util;

namespace hard::tasks::level_set {

// Periodic wrapper
static int
wrapi(int a, int n) {
  int r = a % n;
  return (r < 0) ? (r + n) : r;
}

// Freeze material labels at the start of each step
template<std::size_t D>
void
freeze_material(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> material_a,
  typename field<int>::template accessor<wo, ro> frozen_material_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto mat = m.template mdcolex<is::cells>(material_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      fmat(i) = mat(i);
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(mat,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      fmat(i, j) = mat(i, j);
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(mat,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      fmat(i, j, k) = mat(i, j, k);
    };
  }
}

// upwind HJ RHS (acts on material-signed psi_m)
template<std::size_t D>
void
rhs_upwind(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> velocity_a,
  typename field<int>::template accessor<ro, ro> frozen_material_a,
  typename field<double>::template accessor<wo, ro> rhs_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto vel = m.template mdcolex<is::cells>(velocity_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);

  if constexpr(D == 1) {
    const int nx = m.template size<ax::x, dm::quantities>();

    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const int ip1 = wrapi(i + 1, nx);
      const int im1 = wrapi(i - 1, nx);

      const double up = std::max(0.0, vel(i).x());
      const double um = std::min(0.0, vel(i).x());

      // material-signed neighbor values relative to center frozen material
      const double pm_ip1 = (fmat(ip1) == fmat(i)) ? psi(ip1) : -psi(ip1);
      const double pm_im1 = (fmat(im1) == fmat(i)) ? psi(im1) : -psi(im1);
      const double pm_i = psi(i);

      // adv_x in flux-form upwind; rhs = -adv to match working HJ_upwind_L_2D
      const double adv_x = (pm_i * (up - um) + pm_ip1 * um - pm_im1 * up) /
                           m.template delta<ax::x>();

      rhs(i) = -adv_x;
    };
  }
  else if constexpr(D == 2) {
    const int nx = m.template size<ax::x, dm::quantities>();
    const int ny = m.template size<ax::y, dm::quantities>();

    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      const int ip1 = wrapi(i + 1, nx);
      const int im1 = wrapi(i - 1, nx);
      const int jp1 = wrapi(j + 1, ny);
      const int jm1 = wrapi(j - 1, ny);

      const double upx = std::max(0.0, vel(i, j).x());
      const double umx = std::min(0.0, vel(i, j).x());
      const double upy = std::max(0.0, vel(i, j).y());
      const double umy = std::min(0.0, vel(i, j).y());

      const double pm_ip1 =
        (fmat(ip1, j) == fmat(i, j)) ? psi(ip1, j) : -psi(ip1, j);
      const double pm_im1 =
        (fmat(im1, j) == fmat(i, j)) ? psi(im1, j) : -psi(im1, j);
      const double pm_jp1 =
        (fmat(i, jp1) == fmat(i, j)) ? psi(i, jp1) : -psi(i, jp1);
      const double pm_jm1 =
        (fmat(i, jm1) == fmat(i, j)) ? psi(i, jm1) : -psi(i, jm1);
      const double pm_ij = psi(i, j);

      const double adv_x = (pm_ij * (upx - umx) + pm_ip1 * umx - pm_im1 * upx) /
                           m.template delta<ax::x>();

      const double adv_y = (pm_ij * (upy - umy) + pm_jp1 * umy - pm_jm1 * upy) /
                           m.template delta<ax::y>();

      rhs(i, j) = -(adv_x + adv_y);
    };
  }
  else { // D == 3
    const int nx = m.template size<ax::x, dm::quantities>();
    const int ny = m.template size<ax::y, dm::quantities>();
    const int nz = m.template size<ax::z, dm::quantities>();

    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      const int ip1 = wrapi(i + 1, nx);
      const int im1 = wrapi(i - 1, nx);
      const int jp1 = wrapi(j + 1, ny);
      const int jm1 = wrapi(j - 1, ny);
      const int kp1 = wrapi(k + 1, nz);
      const int km1 = wrapi(k - 1, nz);

      const double upx = std::max(0.0, vel(i, j, k).x());
      const double umx = std::min(0.0, vel(i, j, k).x());
      const double upy = std::max(0.0, vel(i, j, k).y());
      const double umy = std::min(0.0, vel(i, j, k).y());
      const double upz = std::max(0.0, vel(i, j, k).z());
      const double umz = std::min(0.0, vel(i, j, k).z());

      const double pm_ip1 =
        (fmat(ip1, j, k) == fmat(i, j, k)) ? psi(ip1, j, k) : -psi(ip1, j, k);
      const double pm_im1 =
        (fmat(im1, j, k) == fmat(i, j, k)) ? psi(im1, j, k) : -psi(im1, j, k);
      const double pm_jp1 =
        (fmat(i, jp1, k) == fmat(i, j, k)) ? psi(i, jp1, k) : -psi(i, jp1, k);
      const double pm_jm1 =
        (fmat(i, jm1, k) == fmat(i, j, k)) ? psi(i, jm1, k) : -psi(i, jm1, k);
      const double pm_kp1 =
        (fmat(i, j, kp1) == fmat(i, j, k)) ? psi(i, j, kp1) : -psi(i, j, kp1);
      const double pm_km1 =
        (fmat(i, j, km1) == fmat(i, j, k)) ? psi(i, j, km1) : -psi(i, j, km1);
      const double pm_ijk = psi(i, j, k);

      const double adv_x =
        (pm_ijk * (upx - umx) + pm_ip1 * umx - pm_im1 * upx) /
        m.template delta<ax::x>();
      const double adv_y =
        (pm_ijk * (upy - umy) + pm_jp1 * umy - pm_jm1 * upy) /
        m.template delta<ax::y>();
      const double adv_z =
        (pm_ijk * (upz - umz) + pm_kp1 * umz - pm_km1 * upz) /
        m.template delta<ax::z>();

      rhs(i, j, k) = -(adv_x + adv_y + adv_z);
    };
  }
}

// SSPRK3 stage update
template<std::size_t D, std::size_t STEP>
void
advance_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_n0_a, // psi^n
  typename field<double>::template accessor<ro, ro> rhs_a, // k(psi_stage)
  typename field<double>::template accessor<ro, ro>
    psi_prev_a, // stage input (psi1 or psi2)
  typename field<double>::template accessor<wo, ro>
    psi_out_a, // stage output (psi1/psi2/psi3)
  typename single<double>::template accessor<ro> dt_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  static_assert(STEP >= 1 && STEP <= 3, "STEP is 1, 2 or 3");

  auto psi_n0 = m.template mdcolex<is::cells>(psi_n0_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);
  auto psi_prev = m.template mdcolex<is::cells>(psi_prev_a);
  auto psi_out = m.template mdcolex<is::cells>(psi_out_a);

  const double dt = *dt_a;

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if constexpr(STEP == 1) {
        psi_out(i) = psi_n0(i) + dt * rhs(i);
      }
      else if constexpr(STEP == 2) {
        psi_out(i) = 0.75 * psi_n0(i) + 0.25 * (psi_prev(i) + dt * rhs(i));
      }
      else { // STEP == 3
        psi_out(i) =
          (1.0 / 3.0) * psi_n0(i) + (2.0 / 3.0) * (psi_prev(i) + dt * rhs(i));
      }
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      if constexpr(STEP == 1) {
        psi_out(i, j) = psi_n0(i, j) + dt * rhs(i, j);
      }
      else if constexpr(STEP == 2) {
        psi_out(i, j) =
          0.75 * psi_n0(i, j) + 0.25 * (psi_prev(i, j) + dt * rhs(i, j));
      }
      else { // STEP == 3
        psi_out(i, j) = (1.0 / 3.0) * psi_n0(i, j) +
                        (2.0 / 3.0) * (psi_prev(i, j) + dt * rhs(i, j));
      }
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      if constexpr(STEP == 1) {
        psi_out(i, j, k) = psi_n0(i, j, k) + dt * rhs(i, j, k);
      }
      else if constexpr(STEP == 2) {
        psi_out(i, j, k) = 0.75 * psi_n0(i, j, k) +
                           0.25 * (psi_prev(i, j, k) + dt * rhs(i, j, k));
      }
      else { // STEP == 3
        psi_out(i, j, k) =
          (1.0 / 3.0) * psi_n0(i, j, k) +
          (2.0 / 3.0) * (psi_prev(i, j, k) + dt * rhs(i, j, k));
      }
    };
  }
}

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
template<std::size_t D>
void
determine_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> frozen_psi_a, // psi3
  typename field<double>::template accessor<wo, ro> psi_a, // psi^{n+1}
  typename field<int>::template accessor<ro, ro> frozen_material_a, // L^n
  typename field<int>::template accessor<wo, ro> material_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto fpsi = m.template mdcolex<is::cells>(frozen_psi_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);
  auto mat = m.template mdcolex<is::cells>(material_a);

  if constexpr(D == 1) {
    const int nx = m.template size<ax::x, dm::quantities>();

    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(fpsi(i) > 0.0) {
        double psimin = fpsi(i);
        int lmin = fmat(i);

        const int im1 = wrapi(i - 1, nx);
        const int ip1 = wrapi(i + 1, nx);

        if(fmat(im1) != fmat(i) && fpsi(im1) < psimin) {
          psimin = fpsi(im1);
          lmin = fmat(im1);
        }
        if(fmat(ip1) != fmat(i) && fpsi(ip1) < psimin) {
          psimin = fpsi(ip1);
          lmin = fmat(ip1);
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
    const int nx = m.template size<ax::x, dm::quantities>();
    const int ny = m.template size<ax::y, dm::quantities>();

    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      if(fpsi(i, j) > 0.0) {
        double psimin = fpsi(i, j);
        int lmin = fmat(i, j);

        const int ip1 = wrapi(i + 1, nx);
        const int im1 = wrapi(i - 1, nx);
        const int jp1 = wrapi(j + 1, ny);
        const int jm1 = wrapi(j - 1, ny);

        // 4-neighbor scan like working version
        if(fmat(im1, j) != fmat(i, j) && fpsi(im1, j) < psimin) {
          psimin = fpsi(im1, j);
          lmin = fmat(im1, j);
        }
        if(fmat(ip1, j) != fmat(i, j) && fpsi(ip1, j) < psimin) {
          psimin = fpsi(ip1, j);
          lmin = fmat(ip1, j);
        }
        if(fmat(i, jm1) != fmat(i, j) && fpsi(i, jm1) < psimin) {
          psimin = fpsi(i, jm1);
          lmin = fmat(i, jm1);
        }
        if(fmat(i, jp1) != fmat(i, j) && fpsi(i, jp1) < psimin) {
          psimin = fpsi(i, jp1);
          lmin = fmat(i, jp1);
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
    const int nx = m.template size<ax::x, dm::quantities>();
    const int ny = m.template size<ax::y, dm::quantities>();
    const int nz = m.template size<ax::z, dm::quantities>();

    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      if(fpsi(i, j, k) > 0.0) {
        double psimin = fpsi(i, j, k);
        int lmin = fmat(i, j, k);

        const int ip1 = wrapi(i + 1, nx);
        const int im1 = wrapi(i - 1, nx);
        const int jp1 = wrapi(j + 1, ny);
        const int jm1 = wrapi(j - 1, ny);
        const int kp1 = wrapi(k + 1, nz);
        const int km1 = wrapi(k - 1, nz);

        // 6-neighbor scan
        if(fmat(im1, j, k) != fmat(i, j, k) && fpsi(im1, j, k) < psimin) {
          psimin = fpsi(im1, j, k);
          lmin = fmat(im1, j, k);
        }
        if(fmat(ip1, j, k) != fmat(i, j, k) && fpsi(ip1, j, k) < psimin) {
          psimin = fpsi(ip1, j, k);
          lmin = fmat(ip1, j, k);
        }
        if(fmat(i, jm1, k) != fmat(i, j, k) && fpsi(i, jm1, k) < psimin) {
          psimin = fpsi(i, jm1, k);
          lmin = fmat(i, jm1, k);
        }
        if(fmat(i, jp1, k) != fmat(i, j, k) && fpsi(i, jp1, k) < psimin) {
          psimin = fpsi(i, jp1, k);
          lmin = fmat(i, jp1, k);
        }
        if(fmat(i, j, km1) != fmat(i, j, k) && fpsi(i, j, km1) < psimin) {
          psimin = fpsi(i, j, km1);
          lmin = fmat(i, j, km1);
        }
        if(fmat(i, j, kp1) != fmat(i, j, k) && fpsi(i, j, kp1) < psimin) {
          psimin = fpsi(i, j, kp1);
          lmin = fmat(i, j, kp1);
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

} // namespace hard::tasks::level_set
#endif
