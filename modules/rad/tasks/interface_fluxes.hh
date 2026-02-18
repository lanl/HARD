#ifndef HARD_MODULE_RAD_INTERFACE_FLUXES_HH
#define HARD_MODULE_RAD_INTERFACE_FLUXES_HH

#include "../modules/hydro/numerical_algorithms/riemann_solvers.hh"
#include "types.hh"
#include <cstddef>

namespace hard::tasks::rad {

template<std::size_t Dim>
void
compute_interface_fluxes(flecsi::exec::cpu s,
  std::size_t face_axis,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces_vec<Dim>::accessor<ro, ro> uFace_a,
  typename faces<Dim>::accessor<ro, ro> cFace_a,
  typename faces<Dim>::accessor<ro, ro> EradFace_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> EradF_a,
  // time derivative
  field<double>::accessor<wo, ro> dt_radiation_energy_density_a) noexcept {

  auto [uRight, uLeft] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [cRight, cLeft] = faces<Dim>::mdcolex(m, cFace_a);
  auto [EradRight, EradLeft] = faces<Dim>::mdcolex(m, EradFace_a);
  auto EradF = m.template mdcolex<is::cells>(EradF_a);
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);

  using hard::tasks::util::get_mdiota_policy;
  // Compute (1 / dx^i)
  const auto one_over_dx_i = [&m]() {
    if constexpr(Dim == 1) {
      return std::array<double, 1>{1.0 / m.template delta<ax::x>()};
    }
    else if constexpr(Dim == 2) {
      return std::array<double, 2>{
        1.0 / m.template delta<ax::x>(), 1.0 / m.template delta<ax::y>()};
    }
    else {
      return std::array<double, 3>{1.0 / m.template delta<ax::x>(),
        1.0 / m.template delta<ax::y>(),
        1.0 / m.template delta<ax::z>()};
    }
  }();

  if constexpr(Dim == 1) {
#ifdef DEBUG
    if(face_axis != 0) {
      flog_fatal("1D executable can only compute fluxes along X axis");
    }
#endif
    s.executor().forall(i, (m.template cells<ax::x, dm::corrector>())) {

      // min/max characteristic speeds on left
      const double cT{cLeft(i - 1)};
      const double LminT{uLeft(i - 1).x() - cT};
      const double LmaxT{uLeft(i - 1).x() + cT};

      // min/max characteristic speeds on right
      const double cH{cRight(i)};
      const double LminH{uRight(i).x() - cH};
      const double LmaxH{uRight(i).x() + cH};

      const double f_Erad_T{EradLeft(i - 1) * uLeft(i - 1).x()};
      const double f_Erad_H{EradRight(i) * uRight(i).x()};

      // clang-format off
      EradF(i) = numerical_algorithms::advect_conserved<double>(
        EradLeft(i - 1), EradRight(i), f_Erad_T, f_Erad_H,
        LminT, LmaxT, LminH, LmaxH);
      // clang-format on
    }; // forall

    // Store dF^x/dx into du_dt
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_radiation_energy_density(i) +=
        one_over_dx_i[0] * (EradF(i) - EradF(i + 1));
    }; // forall
  }
  else if constexpr(Dim == 2) {
#ifdef DEBUG
    if(face_axis > 1) {
      flog_fatal("2D executable cannot compute fluxes along Z axis");
    }
#endif
    if(face_axis == 0) {
      auto mdpolicy_qc = get_mdiota_policy(EradF,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(ji, mdpolicy_qc) {

        auto [j, i] = ji;

        // min/max characteristic speeds on left
        const double cT = cLeft(i - 1, j);
        const double LminT = uLeft(i - 1, j).x() - cT;
        const double LmaxT = uLeft(i - 1, j).x() + cT;

        // min/max characteristic speeds on right
        const double cH = cRight(i, j);
        const double LminH = uRight(i, j).x() - cH;
        const double LmaxH = uRight(i, j).x() + cH;

        const double f_Erad_T{EradLeft(i - 1, j) * uLeft(i - 1, j).x()};
        const double f_Erad_H{EradRight(i, j) * uRight(i, j).x()};

        // clang-format off
        EradF(i, j) = numerical_algorithms::advect_conserved<double>(
          EradLeft(i - 1, j), EradRight(i, j), f_Erad_T, f_Erad_H, LminT,
            LmaxT, LminH, LmaxH);
        // clang-format on
      }; // forall

      // Store dF^x/dx into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_radiation_energy_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;
        dt_radiation_energy_density(i, j) +=
          one_over_dx_i[0] * (EradF(i, j) - EradF(i + 1, j));
      }; // forall
    }
    else if(face_axis == 1) {
      auto mdpolicy_cq = get_mdiota_policy(EradF,
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_cq) {

        auto [j, i] = ji;

        // min/max characteristic speeds on left
        const double cT = cLeft(i, j - 1);
        const double LminT = uLeft(i, j - 1).y() - cT;
        const double LmaxT = uLeft(i, j - 1).y() + cT;

        // min/max characteristic speeds on right
        const double cH = cRight(i, j);
        const double LminH = uRight(i, j).y() - cH;
        const double LmaxH = uRight(i, j).y() + cH;

        const double f_Erad_T{EradLeft(i, j - 1) * uLeft(i, j - 1).y()};
        const double f_Erad_H{EradRight(i, j) * uRight(i, j).y()};

        // clang-format off
        EradF(i, j) = numerical_algorithms::advect_conserved<double>(
          EradLeft(i, j - 1), EradRight(i, j), f_Erad_T, f_Erad_H,
            LminT, LmaxT, LminH, LmaxH);
        // clang-format on
      }; // forall

      // Store dF^y/dy into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_radiation_energy_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;
        dt_radiation_energy_density(i, j) +=
          one_over_dx_i[1] * (EradF(i, j) - EradF(i, j + 1));

      }; // forall
    }
  }
  else { // Dim == 3

    if(face_axis == 0) {
      const auto mdpolicy_qqc = get_mdiota_policy(EradF,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(kji, mdpolicy_qqc) {

        auto [k, j, i] = kji;

        // min/max characteristic speeds on left
        const double cT = cLeft(i - 1, j, k);
        const double LminT = uLeft(i - 1, j, k).x() - cT;
        const double LmaxT = uLeft(i - 1, j, k).x() + cT;

        // min/max characteristic speeds on right
        const double cH = cRight(i, j, k);
        const double LminH = uRight(i, j, k).x() - cH;
        const double LmaxH = uRight(i, j, k).x() + cH;

        const double f_Erad_T{EradLeft(i - 1, j, k) * uLeft(i - 1, j, k).x()};
        const double f_Erad_H{EradRight(i, j, k) * uRight(i, j, k).x()};

        // clang-format off
        EradF(i, j, k) = numerical_algorithms::advect_conserved<double>(
          EradLeft(i - 1, j, k), EradRight(i, j, k), f_Erad_T, f_Erad_H,
            LminT, LmaxT, LminH, LmaxH);
        // clang-format on

      }; // forall

      // Store dF^x/dx into du_dt
      // Store dF^x/dx into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_radiation_energy_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_radiation_energy_density(i, j, k) +=
          one_over_dx_i[0] * (EradF(i, j, k) - EradF(i + 1, j, k));
      }; // forall
    }
    else if(face_axis == 1) {
      const auto mdpolicy_qcq = get_mdiota_policy(EradF,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qcq) {

        auto [k, j, i] = kji;

        // min/max characteristic speeds on left
        const double cT = cLeft(i, j - 1, k);
        const double LminT = uLeft(i, j - 1, k).y() - cT;
        const double LmaxT = uLeft(i, j - 1, k).y() + cT;

        // min/max characteristic speeds on right
        const double cH = cRight(i, j, k);
        const double LminH = uRight(i, j, k).y() - cH;
        const double LmaxH = uRight(i, j, k).y() + cH;

        const double f_Erad_T{EradLeft(i, j - 1, k) * uLeft(i, j - 1, k).y()};
        const double f_Erad_H{EradRight(i, j, k) * uRight(i, j, k).y()};

        // clang-format off
        EradF(i, j, k) = numerical_algorithms::advect_conserved<double>(
          EradLeft(i, j - 1, k), EradRight(i, j, k), f_Erad_T, f_Erad_H, LminT,
            LmaxT, LminH, LmaxH);
        // clang-format on

      }; // forall

      // Store dF^y/dy into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_radiation_energy_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_radiation_energy_density(i, j, k) +=
          one_over_dx_i[1] * (EradF(i, j, k) - EradF(i, j + 1, k));
      }; // forall
    }
    else {
      const auto mdpolicy_cqq = get_mdiota_policy(EradF,
        m.template cells<ax::z, dm::corrector>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_cqq) {

        auto [k, j, i] = kji;

        // min/max characteristic speeds on left
        const double cT = cLeft(i, j, k - 1);
        const double LminT = uLeft(i, j, k - 1).z() - cT;
        const double LmaxT = uLeft(i, j, k - 1).z() + cT;

        // min/max characteristic speeds on right
        const double cH = cRight(i, j, k);
        const double LminH = uRight(i, j, k).z() - cH;
        const double LmaxH = uRight(i, j, k).z() + cH;

        const double f_Erad_T{EradLeft(i, j, k - 1) * uLeft(i, j, k - 1).z()};
        const double f_Erad_H{EradRight(i, j, k) * uRight(i, j, k).z()};

        // clang-format off
        EradF(i, j, k) = numerical_algorithms::advect_conserved<double>(
          EradLeft(i, j, k - 1), EradRight(i, j, k), f_Erad_T, f_Erad_H,
            LminT, LmaxT, LminH, LmaxH);
        // clang-format on

      }; // forall

      // Store dF^z/dz into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_radiation_energy_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_radiation_energy_density(i, j, k) +=
          one_over_dx_i[2] * (EradF(i, j, k) - EradF(i, j, k + 1));
      }; // forall
    }
  }
}

} // namespace hard::tasks::rad

#endif
