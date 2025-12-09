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
rhs_upwind(const Grid & g,
  const State & stage,
  const int * Lfrozen,
  double t,
  const Vel & vel,
  double * Lout) {
  auto psim = [&](int ic, int jc, int kc, int in, int jn, int kn) {
    int idc = g.idx(ic, jc, kc);
    int idn = g.idx(g.wi(in), g.wj(jn), g.wk(kn));
    return (Lfrozen[idn] == Lfrozen[idc]) ? stage.psi[idn] : -stage.psi[idn];
  };

  for(int k = 0; k < g.nz; ++k)
    for(int j = 0; j < g.ny; ++j)
      for(int i = 0; i < g.nx; ++i) {
        int id = g.idx(i, j, k);
        const double x = (i + 0.0) * g.dx, y = (j + 0.0) * g.dy,
                     z = (k + 0.0) * g.dz;
        double u, v, w;
        vel.sample(x, y, z, t, u, v, w);

        const double pm = stage.psi[id];
        const double im1 = psim(i, j, k, i - 1, j, k),
                     ip1 = psim(i, j, k, i + 1, j, k);
        const double dmx = (pm - im1) / g.dx, dpx = (ip1 - pm) / g.dx;
        double adv = std::max(0.0, u) * dmx + std::min(0.0, u) * dpx;

        if(g.dim >= 2) {
          const double jm1 = psim(i, j, k, i, j - 1, k),
                       jp1 = psim(i, j, k, i, j + 1, k);
          const double dmy = (pm - jm1) / g.dy, dpy = (jp1 - pm) / g.dy;
          double vpos = std::max(0.0, v), vneg = std::min(0.0, v);
          adv += vpos * dmy + vneg * dpy;
        }
        if(g.dim == 3) {
          const double km1 = psim(i, j, k, i, j, k - 1),
                       kp1 = psim(i, j, k, i, j, k + 1);
          const double dmz = (pm - km1) / g.dz, dpz = (kp1 - pm) / g.dz;
          double wpos = std::max(0.0, w), wneg = std::min(0.0, w);
          adv += wpos * dmz + wneg * dpz;
        }

        Lout[id] = -adv; // psi_t = -u · ∇psi
      }
}

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
inline void
determine_psi(const Grid & g,
  const double * psi3,
  const int * Lold,
  double * psi_out,
  int * L_out) {
  if(g.dim == 1) {
    for(int i = 0; i < g.nx; ++i) {
      int id = g.idx(i, 0, 0);
      if(psi3[id] > 0.0) {
        double psimin = psi3[id];
        int lmin = Lold[id];
        int idL = g.idx(g.wi(i - 1), 0, 0), idR = g.idx(g.wi(i + 1), 0, 0);
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

  const int di6[6] = {-1, 1, 0, 0, 0, 0}, dj6[6] = {0, 0, -1, 1, 0, 0},
            dk6[6] = {0, 0, 0, 0, -1, 1};
  const int nnei = (g.dim == 2 ? 4 : 6);

  for(int k = 0; k < g.nz; ++k)
    for(int j = 0; j < g.ny; ++j)
      for(int i = 0; i < g.nx; ++i) {
        int id = g.idx(i, j, k);
        if(psi3[id] <= 0.0) {
          psi_out[id] = psi3[id];
          L_out[id] = Lold[id];
          continue;
        }
        double psimin = psi3[id];
        int lmin = Lold[id];
        for(int n = 0; n < nnei; ++n) {
          int q = g.idx(g.wi(i + di6[n]), g.wj(j + dj6[n]), g.wk(k + dk6[n]));
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
template<class Vel>
inline void
ssprk3_step(const Grid & g, State S, const Vel & vel, double t, double dt) {
  const int N = g.N();
  std::vector<int> Lfrozen(S.label, S.label + N);
  std::vector<double> L(N), psi1(N), psi2(N), psi3(N);

  rhs_upwind(g, S, Lfrozen.data(), t, vel, L.data());
  for(int q = 0; q < N; ++q)
    psi1[q] = S.psi[q] + dt * L[q];

  State S1{psi1.data(), S.label};
  rhs_upwind(g, S1, Lfrozen.data(), t + dt, vel, L.data());
  for(int q = 0; q < N; ++q)
    psi2[q] = 0.75 * S.psi[q] + 0.25 * (psi1[q] + dt * L[q]);

  State S2{psi2.data(), S.label};
  rhs_upwind(g, S2, Lfrozen.data(), t + 0.5 * dt, vel, L.data());
  for(int q = 0; q < N; ++q)
    psi3[q] = (1.0 / 3.0) * S.psi[q] + (2.0 / 3.0) * (psi2[q] + dt * L[q]);

  std::vector<double> psi_next(N);
  std::vector<int> L_next(N);
  determine_psi(g, psi3.data(), S.label, psi_next.data(), L_next.data());
  std::copy(psi_next.begin(), psi_next.end(), S.psi);
  std::copy(L_next.begin(), L_next.end(), S.label);
}

// gradients and unit normals (for GFM & extrapolation)
inline void
gradient_centered(const Grid & g,
  const double * psi,
  double * dpsidx,
  double * dpsidy,
  double * dpsidz) {
  for(int k = 0; k < g.nz; ++k)
    for(int j = 0; j < g.ny; ++j)
      for(int i = 0; i < g.nx; ++i) {
        int id = g.idx(i, j, k);

        auto at = [&](int I, int J, int K) {
          return g.idx(g.wi(I), g.wj(J), g.wk(K));
        };

        int im1 = at(i - 1, j, k), ip1 = at(i + 1, j, k);
        dpsidx[id] = (psi[ip1] - psi[im1]) / (2 * g.dx);

        if(g.dim >= 2) {
          int jm1 = at(i, j - 1, k), jp1 = at(i, j + 1, k);
          dpsidy[id] = (psi[jp1] - psi[jm1]) / (2 * g.dy);
        }
        else
          dpsidy[id] = 0.0;

        if(g.dim == 3) {
          int km1 = at(i, j, k - 1), kp1 = at(i, j, k + 1);
          dpsidz[id] = (psi[kp1] - psi[km1]) / (2 * g.dz);
        }
        else
          dpsidz[id] = 0.0;
      }
}

inline void
unit_normals(const Grid & g,
  const double * dpsidx,
  const double * dpsidy,
  const double * dpsidz,
  double * nx,
  double * ny,
  double * nz,
  double eps = 1e-12) {
  const int N = g.N();
  for(int id = 0; id < N; ++id) {
    double gx = dpsidx[id], gy = dpsidy[id], gz = dpsidz[id];
    double mag = std::sqrt(gx * gx + gy * gy + gz * gz) + eps;
    nx[id] = gx / mag;
    ny[id] = gy / mag;
    nz[id] = gz / mag;
  }
}

// constant-normal extrapolation utility
// Extends scalar q off the interface (one forward Euler step of
// \partial tau q  + sing(\psi) n·\nabla q = 0 ).
template<class SampleQ>
inline void
extend_constant_normal(const Grid & g,
  const double * psi,
  const double * nx,
  const double * ny,
  const double * nz,
  double dtau,
  double * q,
  SampleQ sample) {
  auto at = [&](
              int I, int J, int K) { return g.idx(g.wi(I), g.wj(J), g.wk(K)); };
  std::vector<double> qn(q, q + g.N());

  for(int k = 0; k < g.nz; ++k)
    for(int j = 0; j < g.ny; ++j)
      for(int i = 0; i < g.nx; ++i) {
        int id = g.idx(i, j, k);
        double s = (psi[id] >= 0.0 ? 1.0 : -1.0);
        // first-order upwind along normal: pick upwind neighbor w.r.t. normal
        // components
        int iU = at(i - (nx[id] > 0 ? 1 : 0), j, k),
            iD = at(i + (nx[id] < 0 ? 1 : 0), j, k);
        int jU = at(i, j - (ny[id] > 0 ? 1 : 0), k),
            jD = at(i, j + (ny[id] < 0 ? 1 : 0), k);
        int kU = at(i, j, k - (nz[id] > 0 ? 1 : 0)),
            kD = at(i, j, k + (nz[id] < 0 ? 1 : 0));

        double dqdx = (q[id] - q[iU]) / g.dx;
        double dqdy = (g.dim >= 2 ? (q[id] - q[jU]) / g.dy : 0.0);
        double dqdz = (g.dim == 3 ? (q[id] - q[kU]) / g.dz : 0.0);

        double ndgrad = nx[id] * dqdx + ny[id] * dqdy + nz[id] * dqdz;
        qn[id] = q[id] - dtau * s * ndgrad;
      }
  std::copy(qn.begin(), qn.end(), q);
}

} // namespace spec::level_set::core

#endif
