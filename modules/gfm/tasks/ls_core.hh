#ifndef HARD_MODULES_GFM_TASKS_LS_CORE_HH
#define HARD_MODULES_GFM_TASKS_LS_CORE_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/mesh.hh"
#include "../state.hh"

using namespace hard::tasks::util;

namespace hard::tasks::level_set {

// Freeze material labels at the start of each step
template<std::size_t D>
void
freeze_material(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> material_a,
  typename field<int>::template accessor<wo, ro> frozen_material_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  s.executor().forall(
    i, flecsi::util::iota_view({}, frozen_material_a.span().size())) {
    frozen_material_a(i) = material_a(i);
  };
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
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const int ip1 = i + 1;
      const int im1 = i - 1;

      auto [um, up] = minmax(0.0, vel(i));

      // material-signed neighbor values relative to center frozen material
      const double pm_ip1 = (fmat(ip1) == fmat(i)) ? psi(ip1) : -psi(ip1);
      const double pm_im1 = (fmat(im1) == fmat(i)) ? psi(im1) : -psi(im1);
      const double pm_i = psi(i);

      // adv_x in flux-form upwind; rhs = -adv to match working HJ_upwind_L_2D
      const double adv_x =
        (pm_i * (up.x() - um.x()) + pm_ip1 * um.x() - pm_im1 * up.x()) /
        m.template delta<ax::x>();

      rhs(i) = -adv_x;
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      const int ip1 = i + 1;
      const int im1 = i - 1;
      const int jp1 = j + 1;
      const int jm1 = j - 1;

      auto [um, up] = minmax(0.0, vel(i, j));

      // Compute material-signed neighbors
      std::array<std::tuple<int, int>, 2 * D> neighbor_coords = {
        {{im1, j}, {ip1, j}, {i, jm1}, {i, jp1}}};

      std::array<double, 2 * D> pm;
      for(int n = 0; n < 2 * D; ++n) {
        auto [ni, nj] = neighbor_coords[n];
        pm[n] = (fmat(ni, nj) == fmat(i, j)) ? psi(ni, nj) : -psi(ni, nj);
      }
      const double pm_ij = psi(i, j);

      // pm[0]=im1, pm[1]=ip1, pm[2]=jm1, pm[3]=jp1
      const double adv_x =
        (pm_ij * (up.x() - um.x()) + pm[1] * um.x() - pm[0] * up.x()) /
        m.template delta<ax::x>();

      const double adv_y =
        (pm_ij * (up.y() - um.y()) + pm[3] * um.y() - pm[2] * up.y()) /
        m.template delta<ax::y>();

      rhs(i, j) = -(adv_x + adv_y);
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      const int ip1 = i + 1;
      const int im1 = i - 1;
      const int jp1 = j + 1;
      const int jm1 = j - 1;
      const int kp1 = k + 1;
      const int km1 = k - 1;

      auto [um, up] = minmax(0.0, vel(i, j, k));

      // Compute material-signed neighbors
      std::array<std::tuple<int, int, int>, 2 * D> neighbor_coords = {
        {{im1, j, k},
          {ip1, j, k},
          {i, jm1, k},
          {i, jp1, k},
          {i, j, km1},
          {i, j, kp1}}};

      std::array<double, 2 * D> pm;
      for(int n = 0; n < 2 * D; ++n) {
        auto [ni, nj, nk] = neighbor_coords[n];
        pm[n] = (fmat(ni, nj, nk) == fmat(i, j, k)) ? psi(ni, nj, nk)
                                                    : -psi(ni, nj, nk);
      }
      const double pm_ijk = psi(i, j, k);

      // pm[0]=im1, pm[1]=ip1, pm[2]=jm1, pm[3]=jp1, pm[4]=km1, pm[5]=kp1
      const double adv_x =
        (pm_ijk * (up.x() - um.x()) + pm[1] * um.x() - pm[0] * up.x()) /
        m.template delta<ax::x>();
      const double adv_y =
        (pm_ijk * (up.y() - um.y()) + pm[3] * um.y() - pm[2] * up.y()) /
        m.template delta<ax::y>();
      const double adv_z =
        (pm_ijk * (up.z() - um.z()) + pm[5] * um.z() - pm[4] * up.z()) /
        m.template delta<ax::z>();

      rhs(i, j, k) = -(adv_x + adv_y + adv_z);
    };
  }
}

// SSPRK3 stage update
template<std::size_t D, std::size_t STEP>
void
advance_psi_s1(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_n0_a, // psi^n
  typename field<double>::template accessor<ro, ro> rhs_a, // k(psi_stage)
  typename field<double>::template accessor<wo, ro>
    psi_out_a, // stage output (psi1/psi2/psi3)
  typename single<double>::template accessor<ro> dt_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  s.executor().forall(i, flecsi::util::iota_view({}, psi_out_a.span().size())) {
    psi_out_a(i) = psi_n0_a(i) + *dt_a * rhs_a(i);
  };
}

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

  s.executor().forall(i, flecsi::util::iota_view({}, psi_out_a.span().size())) {
    if constexpr(STEP == 1) {
      psi_out_a(i) = psi_n0_a(i) + *dt_a * rhs_a(i);
    }
    else if constexpr(STEP == 2) {
      psi_out_a(i) =
        0.75 * psi_n0_a(i) + 0.25 * (psi_prev_a(i) + *dt_a * rhs_a(i));
    }
    else { // STEP == 3
      psi_out_a(i) = (1.0 / 3.0) * psi_n0_a(i) +
                     (2.0 / 3.0) * (psi_prev_a(i) + *dt_a * rhs_a(i));
    }
  };
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
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(fpsi(i) > 0.0) {
        double psimin = fpsi(i);
        int lmin = fmat(i);

        // Array of neighbor indices
        const int ii = static_cast<int>(i);
        std::array<int, 2 * D> neighbors = {ii - 1, ii + 1};

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
        const int ii = static_cast<int>(i);
        const int jj = static_cast<int>(j);
        std::array<std::tuple<int, int>, 2 * D> neighbors = {
          {{ii - 1, jj}, {ii + 1, jj}, {ii, jj - 1}, {ii, jj + 1}}};

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
        const int ii = static_cast<int>(i);
        const int jj = static_cast<int>(j);
        const int kk = static_cast<int>(k);
        std::array<std::tuple<int, int, int>, 2 * D> neighbors = {
          {{ii - 1, jj, kk},
            {ii + 1, jj, kk},
            {ii, jj - 1, kk},
            {ii, jj + 1, kk},
            {ii, jj, kk - 1},
            {ii, jj, kk + 1}}};

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

} // namespace hard::tasks::level_set
#endif
