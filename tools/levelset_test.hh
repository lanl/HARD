// TODO: part of unit test?

#ifndef LEVELSET_TEST_HH
#define LEVELSET_TEST_HH

#include "../labels.hh"
#include "../state.hh"
#include "../types.hh"

#include "../state.hh"
#include "spec/levelset/ls_core.hh"
#include "spec/levelset/ls_initial_velocity_test.hh"
#include "spec/mesh.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace hard::tasks::level_set_test {

template<std::size_t D, class Scheduler>
void
levelset_time_reversal(flecsi::exec::accelerator s,
  Scheduler & sched,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<rw> psi_a,
  typename field<vec<D>>::template accessor<ro> velocity_a,
  typename field<int>::template accessor<rw> label_a,
  typename mesh<D>::topology & topo, // TODO: only for standalone IO
  double T,
  double CFL,
  int num_steps) {

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto u = m.template mdcolex<is::cells>(velocity_a);
  auto label = m.template mdcolex<is::cells>(label_a);

  // --- grid meta (host-side scalars) ---
  int nx = 1, ny = 1, nz = 1;
  double dx = 1.0, dy = 1.0, dz = 1.0;

  if constexpr(D == 1) {
    nx = m.template size<ax::x, dm::quantities>();
    ny = 1;
    nz = 1;
    dx = m.template delta<ax::x>();
    dy = 1.0;
    dz = 1.0;
  }
  else if constexpr(D == 2) {
    nx = m.template size<ax::x, dm::quantities>();
    ny = m.template size<ax::y, dm::quantities>();
    nz = 1;
    dx = m.template delta<ax::x>();
    dy = m.template delta<ax::y>();
    dz = 1.0;
  }
  else /* D == 3 */ {
    nx = m.template size<ax::x, dm::quantities>();
    ny = m.template size<ax::y, dm::quantities>();
    nz = m.template size<ax::z, dm::quantities>();
    dx = m.template delta<ax::x>();
    dy = m.template delta<ax::y>();
    dz = m.template delta<ax::z>();
  }

  // --- interface geometry : simple circular ---
  const double Lx = dx * nx;
  const double Ly = dy * ny;
  const double Lz = dz * nz;

  const double cx = 0.5 * Lx;
  const double cy = (D >= 2 ? 0.5 * Ly : 0.0);
  const double cz = (D == 3 ? 0.5 * Lz : 0.0);

  const double R0 = 0.2 * Lx;

  // --- Initialization: set psi + label ---
  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double x = m.template center<ax::x>(i);
      const double r = std::abs(x - cx);

      psi(i) = r - R0;
      label(i) = (psi(i) <= 0.0 ? 1 : 0);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      const double x = m.template center<ax::x>(i);
      const double y = m.template center<ax::y>(j);
      const double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));

      psi(i, j) = r - R0;
      label(i, j) = (psi(i, j) <= 0.0 ? 1 : 0);
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      const double x = m.template center<ax::x>(i);
      const double y = m.template center<ax::y>(j);
      const double z = m.template center<ax::z>(k);

      const double r = std::sqrt(
        (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));

      psi(i, j, k) = r - R0;
      label(i, j, k) = (psi(i, j, k) <= 0.0 ? 1 : 0);
    };
  }

  // --- time loop (0 -> 2T) ---
  double hmin = dx;
  if constexpr(D >= 2)
    hmin = std::min(hmin, dy);
  if constexpr(D == 3)
    hmin = std::min(hmin, dz);

  // simple dt TODO: tighten using max|u| from velocity field with hydro
  // integration
  double dt = CFL * hmin;

  double t = 0.0;
  const double t_end = 2.0 * T;

  int step = 0;
  while(t < t_end && step < num_steps) {
    if(t + dt > t_end)
      dt = t_end - t;

    // Fill velocity_a at time t
    // In this test, we use time reversal test as described in Aslam et al. 2025
    // TODO : test with hydro velocity

    auto vel = m.template mdcolex<is::cells>(velocity_a);
    if constexpr(D == 1) {
      spec::level_set::vel::TimeReversal1D vfoo{T};

      s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
        const double x = m.template center<ax::x>(i);

        double ux = 0.0, uy = 0.0, uz = 0.0;
        vfoo.sample(x, 0.0, 0.0, t, ux, uy, uz);

        vel(i).x() = ux;
      }; // forall
    }
    else if constexpr(D == 2) {
      spec::level_set::vel::TimeReversal2D vfoo{T};

      auto mdpolicy_qq = get_mdiota_policy(vel,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;

        const double x = m.template center<ax::x>(i);
        const double y = m.template center<ax::y>(j);

        double ux = 0.0, uy = 0.0, uz = 0.0;
        vfoo.sample(x, y, 0.0, t, ux, uy, uz);

        vel(i, j).x() = ux;
        vel(i, j).y() = uy;
      }; // forall
    }
    else /* D == 3 */ {
      spec::level_set::vel::TimeReversal3D vfoo{T};

      auto mdpolicy_qqq = get_mdiota_policy(vel,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;

        const double x = m.template center<ax::x>(i);
        const double y = m.template center<ax::y>(j);
        const double z = m.template center<ax::z>(k);

        double ux = 0.0, uy = 0.0, uz = 0.0;
        vfoo.sample(x, y, z, t, ux, uy, uz);

        vel(i, j, k).x() = ux;
        vel(i, j, k).y() = uy;
        vel(i, j, k).z() = uz;
      }; // forall
    }

    hard::tasks::level_set::StepParams P{t, dt};

    sched.execute(hard::task::level_set::advance_levelset_task<D>,
      m,
      psi_a,
      label_a,
      velocity_a,
      P);

    t += dt;
    ++step;

    if(step % 10 == 0) {
      hard::task::level_set::dump_levelset_snapshot<mesh<D>>(
        topo, "ls_tr", step);
    }
  }
}

} // namespace hard::tasks::level_set_test

#endif
