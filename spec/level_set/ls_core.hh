#ifndef SPEC_LS_CORE_HH
#define SPEC_LS_CORE_HH

#include "../../app/state.hh"
#include "../mesh.hh"
#include "app/types.hh"
#include "spec/labels.hh"

namespace spec::level_set::core {

namespace task {

// upwind HJ RHS (acts on material-signed psi_m)
template<std::size_t D>
void
rhs_upwind(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> psi_a,
  typename field<vec<D>>::template accessor<ro> velocity_a,
  typename field<int>::template accessor<rw> frozen_label_a,
  typename field<int>::template accessor<rw> label_a) {

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto u = m.template mdcolex<is::cells>(velocity_a);
  auto label = m.template mdcolex<is::cells>(label_a);
  auto f_label = m.template mdcolex<is::cells>(frozen_label_a);

  // TODO: The velocity is not sampled right now at a timestep
  // This is why the dt was being added before, but I removed it from
  // the signature for now

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double max_u_x = std::max(0.0, u(i).x());
      const double min_u_x = std::min(0.0, u(i).x());
      const double psi_ip1 =
        (f_label(i) == f_label(i + 1)) ? psi(i + 1) : -psi(i + 1);
      const double psi_im1 =
        (f_label(i) == f_label(i - 1)) ? psi(i - 1) : -psi(i - 1);

      label(i) -=
        (psi(i) * (max_u_x - min_u_x) + psi_ip1 * min_u_x - psi_im1 * max_u_x) /
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
      const double psi_ip1 =
        (f_label(i, j) == f_label(i + 1, j)) ? psi(i + 1, j) : -psi(i + 1, j);
      const double psi_im1 =
        (f_label(i, j) == f_label(i - 1, j)) ? psi(i - 1, j) : -psi(i - 1, j);
      const double psi_jp1 =
        (f_label(i, j) == f_label(i, j + 1)) ? psi(i, j + 1) : -psi(i, j + 1);
      const double psi_jm1 =
        (f_label(i, j) == f_label(i, j - 1)) ? psi(i, j - 1) : -psi(i, j - 1);

      label(i, j) -= (psi(i, j) * (max_u_x - min_u_x) + psi_ip1 * min_u_x -
                       psi_im1 * max_u_x) /
                       m.template delta<ax::x>() +
                     (psi(i, j) * (max_u_y - min_u_y) + psi_jp1 * min_u_y -
                       psi_jm1 * max_u_y) /
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
      const double psi_ip1 = (f_label(i, j, k) == f_label(i + 1, j, k))
                               ? psi(i + 1, j, k)
                               : -psi(i + 1, j, k);
      const double psi_im1 = (f_label(i, j, k) == f_label(i - 1, j, k))
                               ? psi(i - 1, j, k)
                               : -psi(i - 1, j, k);
      const double psi_jp1 = (f_label(i, j, k) == f_label(i, j + 1, k))
                               ? psi(i, j + 1, k)
                               : -psi(i, j + 1, k);
      const double psi_jm1 = (f_label(i, j, k) == f_label(i, j - 1, k))
                               ? psi(i, j - 1, k)
                               : -psi(i, j - 1, k);
      const double psi_kp1 = (f_label(i, j, k) == f_label(i, j, k + 1))
                               ? psi(i, j, k + 1)
                               : -psi(i, j, k + 1);
      const double psi_km1 = (f_label(i, j, k) == f_label(i, j, k - 1))
                               ? psi(i, j, k - 1)
                               : -psi(i, j, k - 1);

      label(i, j, k) -= (psi(i, j, k) * (max_u_x - min_u_x) +
                          psi_ip1 * min_u_x - psi_im1 * max_u_x) /
                          m.template delta<ax::x>() +
                        (psi(i, j, k) * (max_u_y - min_u_y) +
                          psi_jp1 * min_u_y - psi_jm1 * max_u_y) /
                          m.template delta<ax::y>() +
                        (psi(i, j, k) * (max_u_z - min_u_z) +
                          psi_kp1 * min_u_z - psi_km1 * max_u_z) /
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
  typename field<int>::template accessor<rw> frozen_label_a,
  typename field<int>::template accessor<wo> label_a) {

  auto f_psi = m.template mdcolex<is::cells>(frozen_psi_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto f_label = m.template mdcolex<is::cells>(frozen_label_a);
  auto label = m.template mdcolex<is::cells>(label_a);

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
        }
        psi(i, j, k) = (psimin >= 0.0 ? 1.0 : -1.0) * f_psi(i, j, k);
        label(i, j, k) = lmin;
      }
      else {
        psi(i, j, k) = f_psi(i, j, k);
        label(i, j, k) = f_label(i, j, k);
      }
    }; // forall
  }
} // determine_psi

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
template<std::size_t D, std::size_t STEP>
void
advance_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro> psi_a,
  typename field<int>::template accessor<ro> label_a,
  typename field<double>::template accessor<ro> psi_prev_a,
  typename field<double>::template accessor<wo> psi_n_a,
  const double dt) {

  static_assert(STEP >= 1 && STEP <= 3, "STEP is 1, 2 or 3");

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto label = m.template mdcolex<is::cells>(label_a);
  auto psi_prev = m.template mdcolex<is::cells>(psi_prev_a);
  auto psi_n = m.template mdcolex<is::cells>(psi_n_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if constexpr(STEP == 1) {
        psi_n(i) = psi(i) + dt * label(i);
      }
      else if constexpr(STEP == 2) {
        psi_n(i) = (3 * psi(i) + (psi_prev(i) + dt * label(i))) * 0.25;
      }
      else /* STEP == 3 */ {
        psi_n(i) = (psi(i) + 2 * (psi_prev(i) + dt * label(i))) / 3.0;
      }
    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      if constexpr(STEP == 1) {
        psi_n(i, j) = psi(i, j) + dt * label(i, j);
      }
      else if constexpr(STEP == 2) {
        psi_n(i, j) =
          (3 * psi(i, j) + (psi_prev(i, j) + dt * label(i, j))) * 0.25;
      }
      else /* STEP == 3 */ {
        psi_n(i, j) =
          (psi(i, j) + 2 * (psi_prev(i, j) + dt * label(i, j))) / 3.0;
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
      if constexpr(STEP == 1) {
        psi_n(i, j, k) = psi(i, j, k) + dt * label(i, j, k);
      }
      else if constexpr(STEP == 2) {
        psi_n(i, j, k) =
          (3 * psi(i, j, k) + (psi_prev(i, j, k) + dt * label(i, j, k))) * 0.25;
      }
      else /* STEP == 3 */ {
        psi_n(i, j, k) =
          (psi(i, j, k) + 2 * (psi_prev(i, j, k) + dt * label(i, j, k))) / 3.0;
      }
    }; // forall
  }
} // advance_psi
} // namespace task

// SSPRK(3,3) driver with frozen labels
template<std::size_t D>
void
ssprk3_step(control_policy<hard::state, D> & cp) {
  auto & s = cp.state();

  flecsi::scheduler & sc = cp.scheduler();

  sc.execute<task::rhs_upwind<D>>(flecsi::exec::on,
    *s.m,
    s.psi(*s.m),
    s.velocity(*s.m),
    s.frozen_label(*s.m),
    s.label(*s.m));
  sc.execute<task::advance_psi<D, 1>>(flecsi::exec::on,
    *s.m,
    s.psi(*s.m),
    s.label(*s.m),
    s.psi2(*s.m),
    s.psi1(*s.m),
    s.dt(s.gt));

  sc.execute<task::rhs_upwind<D>>(flecsi::exec::on,
    *s.m,
    s.psi1(*s.m),
    s.velocity(*s.m),
    s.frozen_label(*s.m),
    s.label(*s.m));
  sc.execute<task::advance_psi<D, 2>>(flecsi::exec::on,
    *s.m,
    s.psi(*s.m),
    s.label(*s.m),
    s.psi1(*s.m),
    s.psi2(*s.m),
    s.dt(s.gt));

  sc.execute<task::rhs_upwind<D>>(flecsi::exec::on,
    *s.m,
    s.psi2(*s.m),
    s.velocity(*s.m),
    s.frozen_label(*s.m),
    s.label(*s.m));
  sc.execute<task::advance_psi<D, 3>>(flecsi::exec::on,
    *s.m,
    s.psi(*s.m),
    s.label(*s.m),
    s.psi2(*s.m),
    s.psi3(*s.m),
    s.dt(s.gt));

  sc.execute<task::determine_psi<D>>(flecsi::exec::on,
    *s.m,
    s.psi3(*s.m),
    s.psi(*s.m),
    s.frozen_label(*s.m),
    s.label(*s.m));
}
} // namespace spec::level_set::core

#endif
