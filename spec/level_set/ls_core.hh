#ifndef SPEC_LS_CORE_HH
#define SPEC_LS_CORE_HH

#include "../mesh.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace spec::level_set::core {

// Helpers to get sizes and deltas from spec/mesh.hh
template<std::size_t D, class MeshAcc>
inline void
mesh_metrics(const MeshAcc & m,
  int & nx,
  int & ny,
  int & nz,
  double & dx,
  double & dy,
  double & dz) {
  using ax::x;
  using ax::y;
  using ax::z;
  using dm::quantities;
  if constexpr(D == 1) {
    nx = m.template size<x, quantities>();
    ny = 1;
    nz = 1;
    dx = m.template delta<x>();
    dy = 1.0;
    dz = 1.0;
  }
  else if constexpr(D == 2) {
    nx = m.template size<x, quantities>();
    ny = m.template size<y, quantities>();
    nz = 1;
    dx = m.template delta<x>();
    dy = m.template delta<y>();
    dz = 1.0;
  }
  else { // D==3
    nx = m.template size<x, quantities>();
    ny = m.template size<y, quantities>();
    nz = m.template size<z, quantities>();
    dx = m.template delta<x>();
    dy = m.template delta<y>();
    dz = m.template delta<z>();
  }
}

// upwind HJ RHS (acts on material-signed psi_m)
template<std::size_t D, class MeshAcc, class PsiCont, class Vel>
inline void
rhs_upwind(const MeshAcc & m,
  const PsiCont & psi,
  const std::vector<int> & Lfrozen,
  double t,
  const Vel & vel,
  std::vector<double> & Lout) {

  int nx, ny, nz;
  double dx, dy, dz;
  mesh_metrics<D>(m, nx, ny, nz, dx, dy, dz);
  const int N = nx * ny * nz;
  if((int)Lout.size() != N)
    Lout.assign(N, 0.0);

  auto idx = [=](int i, int j, int k) {
    auto wrap = [](int I, int n) {
      I %= n;
      if(I < 0)
        I += n;
      return I;
    };
    i = wrap(i, nx);
    j = wrap(j, ny);
    k = wrap(k, nz);
    return i + nx * (j + ny * k);
  };

  auto psim = [&](int ic, int jc, int kc, int in, int jn, int kn) {
    int idc = idx(ic, jc, kc);
    int idn = idx(in, jn, kn);
    return (Lfrozen[idn] == Lfrozen[idc]) ? stage.psi[idn] : -stage.psi[idn];
  };

  for(int k = 0; k < nz; ++k)
    for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i) {
        int id = idx(i, j, k);

        // map (i,j,k) to physical x,y,z with mesh centers:
        // we assume that mesh is unifrom::
        // x ~ (i+0.5)*dx;
        double x = (i + 0.5) * dx;
        double y = (D >= 2 ? (j + 0.5) * dy : 0.0);
        double z = (D == 3 ? (k + 0.5) * dz : 0.0);

        double u, v, w;
        vel.sample(x, y, z, t, u, v, w);

        const double pm = psi[id];
        const double im1 = psim(i, j, k, i - 1, j, k),
                     ip1 = psim(i, j, k, i + 1, j, k);
        const double dmx = (pm - im1) / dx, dpx = (ip1 - pm) / dx;
        double adv = std::max(0.0, u) * dmx + std::min(0.0, u) * dpx;

        if constexpr(D >= 2) {
          const double jm1 = psim(i, j, k, i, j - 1, k),
                       jp1 = psim(i, j, k, i, j + 1, k);
          const double dmy = (pm - jm1) / dy, dpy = (jp1 - pm) / dy;
          double vpos = std::max(0.0, v), vneg = std::min(0.0, v);
          adv += vpos * dmy + vneg * dpy;
        }
        if constexpr(D == 3) {
          const double km1 = psim(i, j, k, i, j, k - 1),
                       kp1 = psim(i, j, k, i, j, k + 1);
          const double dmz = (pm - km1) / dz, dpz = (kp1 - pm) / dz;
          double wpos = std::max(0.0, w), wneg = std::min(0.0, w);
          adv += wpos * dmz + wneg * dpz;
        }

        Lout[id] = -adv; // psi_t = -u · grad_psi
      }
}

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
template<std::size_t D,
  class Psi3Cont,
  class LabelCont,
  class PsiOutCont,
  class LabelOutCont>
inline void
determine_psi(int nx,
  int ny,
  int nz,
  const Psi3Cont & psi3,
  const LabelCont & Lold,
  PsiOutCont & psi_out,
  LabelOutCont & L_out) {

  auto idx = [=](int i, int j, int k) {
    auto wrap = [](int I, int n) {
      I %= n;
      if(I < 0)
        I += n;
      return I;
    };
    i = wrap(i, nx);
    j = wrap(j, ny);
    k = wrap(k, nz);
    return i + nx * (j + ny * k);
  };

  if constexpr(D == 1) {
    for(int i = 0; i < nx; ++i) {
      int id = idx(i, 0, 0);
      if(psi3[id] > 0.0) {
        double psimin = psi3[id];
        int lmin = Lold[id];
        int idL = idx(i - 1, 0, 0), idR = idx(i + 1, 0, 0);
        if(Lold[idL] != Lold[id] && psi3[idL] < psimin) {
          psimin = psi3[idL];
          lmin = Lold[idL];
        }
        if(Lold[idR] != Lold[id] && psi3[idR] < psimin) {
          psimin = psi3[idR];
          lmin = Lold[idR];
        }
        psi_out[id] = (psimin >= 0.0 ? 1.0 : -1.0) * psi3[id];
        L_out[id] = lmin;
      }
      else {
        psi_out[id] = psi3[id];
        L_out[id] = Lold[id];
      }
    }
    return;
  }

  const int di[6] = {-1, 1, 0, 0, 0, 0}, dj[6] = {0, 0, -1, 1, 0, 0},
            dk[6] = {0, 0, 0, 0, -1, 1};
  const int nnei = (D == 2 ? 4 : 6);

  for(int k = 0; k < nz; ++k)
    for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i) {
        int id = idx(i, j, k);
        if(psi3[id] <= 0.0) {
          psi_out[id] = psi3[id];
          L_out[id] = Lold[id];
          continue;
        }
        double psimin = psi3[id];
        int lmin = Lold[id];
        for(int n = 0; n < nnei; ++n) {
          int q = idx(i + di[n], j + dj[n], k + dk[n]);
          if(Lold[q] != Lold[id] && psi3[q] < psimin) {
            psimin = psi3[q];
            lmin = Lold[q];
          }
        }
        psi_out[id] = (psimin >= 0.0 ? 1.0 : -1.0) * psi3[id];
        L_out[id] = lmin;
      }
}

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
