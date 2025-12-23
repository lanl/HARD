#ifndef SPEC_LS_CORE_HH
#define SPEC_LS_CORE_HH

#include "../mesh.hh"
#include "app/types.hh"
#include "spec/labels.hh"

namespace spec::level_set::core {

namespace hard::tasks {

// upwind HJ RHS (acts on material-signed psi_m)
template<std::size_t D>
void
rhs_upwind(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> psi_a,
  typename field<vec<D>>::template accessor<ro> velocity_a,
  typename field<int>::template accessor<rw> label_a) {

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto u = m.template mdcolex<is::cells>(velocity_a);
  auto label = m.template mdcolex<is::cells>(label_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double max_u_x = std::max(0.0, u(i).x());
      const double min_u_x = std::min(0.0, u(i).x());
      label(i) -= (psi(i) * (max_u_x - min_u_x) + psi(i + 1) * min_u_x -
                    psi(i - 1) * max_u_x) /
                  m.template delta<ax::x>();
    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      const double max_u_x = std::max(0.0, u(i, j).x());
      const double min_u_x = std::min(0.0, u(i, j).x());
      const double max_u_y = std::max(0.0, u(i, j).y());
      const double min_u_y = std::min(0.0, u(i, j).y());
      label(i, j) -= (psi(i, j) * (max_u_x - min_u_x) +
                       psi(i + 1, j) * min_u_x - psi(i - 1, j) * max_u_x) /
                       m.template delta<ax::x>() +
                     (psi(i, j) * (max_u_y - min_u_y) +
                       psi(i, j + 1) * min_u_y - psi(i, j - 1) * max_u_y) /
                       m.template delta<ax::y>();
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      const double max_u_x = std::max(0.0, u(i, j, k).x());
      const double min_u_x = std::min(0.0, u(i, j, k).x());
      const double max_u_y = std::max(0.0, u(i, j, k).y());
      const double min_u_y = std::min(0.0, u(i, j, k).y());
      const double max_u_z = std::max(0.0, u(i, j, k).z());
      const double min_u_z = std::min(0.0, u(i, j, k).z());
      label(i, j, k) -=
        (psi(i, j, k) * (max_u_x - min_u_x) + psi(i + 1, j, k) * min_u_x -
          psi(i - 1, j, k) * max_u_x) /
          m.template delta<ax::x>() +
        (psi(i, j, k) * (max_u_y - min_u_y) + psi(i, j + 1, k) * min_u_y -
          psi(i, j - 1, k) * max_u_y) /
          m.template delta<ax::y>() +
        (psi(i, j, k) * (max_u_z - min_u_z) + psi(i, j, k + 1) * min_u_z -
          psi(i, j, k - 1) * max_u_z) /
          m.template delta<ax::z>();
    }; // forall
  }
}

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
template<std::size_t D>
void
determine_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> frozen_psi_a,
  typename field<double>::template accessor<wo> psi_a,
  typename field<int>::template accessor<ro> frozen_label_a,
  typename field<int>::template accessor<wo> label_a) {

  auto f_psi = m.template mdcolex<is::cells>(frozen_psi_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto f_label = m.template mdcolex<is::cells>(frozen_label_a);
  auto Label = m.template mdcolex<is::cells>(label_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(f_psi(i) > 0) {
        double psimin = f_psi(i);
        int lmin = f_label(i);
        for(int pm = -1; pm <= 1; pm += 2) {
          if(f_label(i + pm) != f_label(i) && f_psi(i + pm) < psimin) {
            psimin = f_psi(i + pm);
            lmin = f_label(i + pm);
          }
        }
        psi(i) = (psimin >= 0.0 ? 1.0 : -1.0) * f_psi(i);
        label(i) = lmin;
      }
      else {
        psi(i) = f_psi(i);
        label(i) = f_label(i);
      }
    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      if(f_psi(i, j) > 0) {
        double psimin = f_psi(i, j);
        int lmin = f_label(i, j);
        for(int pm = -1; pm <= 1; pm += 2) {
          if(f_label(i + pm, j) != f_label(i, j) && f_psi(i + pm, j) < psimin) {
            psimin = f_psi(i + pm, j);
            lmin = f_label(i + pm, j);
          }
          if(f_label(i, j + pm) != f_label(i, j) && f_psi(i, j + pm) < psimin) {
            psimin = f_psi(i, j + pm);
            lmin = f_label(i, j + pm);
          }
        }
        psi(i, j) = (psimin >= 0.0 ? 1.0 : -1.0) * f_psi(i, j);
        label(i, j) = lmin;
      }
      else {
        psi(i, j) = f_psi(i, j);
        label(i, j) = f_label(i, j);
      }
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      if(f_psi(i, j, k) > 0) {
        double psimin = f_psi(i, j, k);
        int lmin = f_label(i, j, k);
        for(int pm = -1; pm <= 1; pm += 2) {
          if(f_label(i + pm, j, k) != f_label(i, j, k) &&
             f_psi(i + pm, j, k) < psimin) {
            psimin = f_psi(i + pm, j, k);
            lmin = f_label(i + pm, j, k);
          }
          if(f_label(i, j + pm, k) != f_label(i, j, k) &&
             f_psi(i, j + pm, k) < psimin) {
            psimin = f_psi(i, j + pm, k);
            lmin = f_label(i, j + pm, k);
          }
          if(f_label(i, j, k + pm) != f_label(i, j, k) &&
             f_psi(i, j, k + pm) < psimin) {
            psimin = f_psi(i, j, k + pm);
            lmin = f_label(i, j, k + pm);
          }
          psi(i, j, k) = (psimin >= 0.0 ? 1.0 : -1.0) * f_psi(i, j, k);
          label(i, j, k) = lmin;
        }
      }
      else {
        psi(i, j, k) = f_psi(i, j, k);
        label(i, j, k) = f_label(i, j, k);
      }
    }; // forall
  }

} // namespace hard::tasks

// SSPRK(3,3) driver with frozen labels
template<std::size_t D, class MeshAcc, class PsiAcc, class LabelAcc, class Vel>
inline void
ssprk3_step(const MeshAcc & m,
  PsiAcc & psi_a,
  LabelAcc & lab_a,
  const Vel & vel,
  double t,
  double dt) {

  int nx, ny, nz;
  double dx, dy, dz;
  mesh_metrics<D>(m, nx, ny, nz, dx, dy, dz);
  const int N = nx * ny * nz;

  auto psi_span = psi_a.span();
  auto label_span = lab_a.span();

  std::vector<double> psi(psi_span.begin(), psi_span.end());
  std::vector<int> label(label_span.begin(), label_span.end());

  std::vector<int> Lfrozen = label;
  std::vector<double> L(N), psi1(N), psi2(N), psi3(N);

  rhs_upwind<D>(m, psi, Lfrozen, t, vel, L);
  for(int q = 0; q < N; ++q)
    psi1[q] = psi[q] + dt * L[q];

  rhs_upwind<D>(m, psi1, Lfrozen, t + dt, vel, L);
  for(int q = 0; q < N; ++q)
    psi2[q] = 0.75 * psi[q] + 0.25 * (psi1[q] + dt * L[q]);

  rhs_upwind<D>(m, psi2, Lfrozen, t + 0.5 * dt, vel, L);
  for(int q = 0; q < N; ++q)
    psi3[q] = (1.0 / 3.0) * psi[q] + (2.0 / 3.0) * (psi2[q] + dt * L[q]);

  std::vector<double> psi_next(N);
  std::vector<int> L_next(N);
  determine_psi<D>(nx, ny, nz, psi3, label, psi_next, L_next);

  for(int q = 0; q < N; ++q) {
    psi_span[q] = psi_next[q];
    label_span[q] = L_next[q];
  }
}

} // namespace spec::level_set::core

#endif
