#ifndef HARD_MODULE_RAD_INTERFACE_FLUXES_HH
#define HARD_MODULE_RAD_INTERFACE_FLUXES_HH

#include "../modules/hydro/numerical_algorithms/riemann_solvers.hh"
#include "types.hh"

namespace hard::tasks::rad {

template<std::size_t Dim>
void
compute_interface_fluxes(flecsi::exec::accelerator s,
  std::size_t face_axis,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces_vec<Dim>::template accessor<ro, ro> uFace_a,
  typename faces<Dim>::template accessor<ro, ro> cFace_a,
  typename faces<Dim>::template accessor<ro, ro> erad_face_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> erad_f_a,
  // time derivative
  field<double>::accessor<wo, ro> dt_radiation_energy_density_a) noexcept {

  auto [u_right, u_left] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [c_right, c_left] = faces<Dim>::mdcolex(m, cFace_a);
  auto [erad_right, erad_left] = faces<Dim>::mdcolex(m, erad_face_a);
  auto erad_f = m.template mdcolex<is::cells>(erad_f_a);
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
      const double cT{c_left(i - 1)};
      const double LminT{u_left(i - 1).x() - cT};
      const double LmaxT{u_left(i - 1).x() + cT};

      // min/max characteristic speeds on right
      const double cH{c_right(i)};
      const double LminH{u_right(i).x() - cH};
      const double LmaxH{u_right(i).x() + cH};

      const double f_Erad_T{erad_left(i - 1) * u_left(i - 1).x()};
      const double f_Erad_H{erad_right(i) * u_right(i).x()};

      // clang-format off
      erad_f(i) = numerical_algorithms::advect_conserved<double>(
        erad_left(i - 1), erad_right(i), f_Erad_T, f_Erad_H,
        LminT, LmaxT, LminH, LmaxH);
      // clang-format on
    }; // forall

    // Store dF^x/dx into du_dt
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_radiation_energy_density(i) +=
        one_over_dx_i[0] * (erad_f(i) - erad_f(i + 1));
    }; // forall
  }
  else if constexpr(Dim == 2) {
#ifdef DEBUG
    if(face_axis > 1) {
      flog_fatal("2D executable cannot compute fluxes along Z axis");
    }
#endif
    if(face_axis == 0) {
      auto mdpolicy_qc = get_mdiota_policy(erad_f,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(ji, mdpolicy_qc) {

        auto [j, i] = ji;

        // min/max characteristic speeds on left
        const double cT = c_left(i - 1, j);
        const double LminT = u_left(i - 1, j).x() - cT;
        const double LmaxT = u_left(i - 1, j).x() + cT;

        // min/max characteristic speeds on right
        const double cH = c_right(i, j);
        const double LminH = u_right(i, j).x() - cH;
        const double LmaxH = u_right(i, j).x() + cH;

        const double f_Erad_T{erad_left(i - 1, j) * u_left(i - 1, j).x()};
        const double f_Erad_H{erad_right(i, j) * u_right(i, j).x()};

        // clang-format off
        erad_f(i, j) = numerical_algorithms::advect_conserved<double>(
          erad_left(i - 1, j), erad_right(i, j), f_Erad_T, f_Erad_H, LminT,
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
          one_over_dx_i[0] * (erad_f(i, j) - erad_f(i + 1, j));
      }; // forall
    }
    else if(face_axis == 1) {
      auto mdpolicy_cq = get_mdiota_policy(erad_f,
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_cq) {

        auto [j, i] = ji;

        // min/max characteristic speeds on left
        const double cT = c_left(i, j - 1);
        const double LminT = u_left(i, j - 1).y() - cT;
        const double LmaxT = u_left(i, j - 1).y() + cT;

        // min/max characteristic speeds on right
        const double cH = c_right(i, j);
        const double LminH = u_right(i, j).y() - cH;
        const double LmaxH = u_right(i, j).y() + cH;

        const double f_Erad_T{erad_left(i, j - 1) * u_left(i, j - 1).y()};
        const double f_Erad_H{erad_right(i, j) * u_right(i, j).y()};

        // clang-format off
        erad_f(i, j) = numerical_algorithms::advect_conserved<double>(
          erad_left(i, j - 1), erad_right(i, j), f_Erad_T, f_Erad_H,
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
          one_over_dx_i[1] * (erad_f(i, j) - erad_f(i, j + 1));

      }; // forall
    }
  }
  else { // Dim == 3

    if(face_axis == 0) {
      const auto mdpolicy_qqc = get_mdiota_policy(erad_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(kji, mdpolicy_qqc) {

        auto [k, j, i] = kji;

        // min/max characteristic speeds on left
        const double cT = c_left(i - 1, j, k);
        const double LminT = u_left(i - 1, j, k).x() - cT;
        const double LmaxT = u_left(i - 1, j, k).x() + cT;

        // min/max characteristic speeds on right
        const double cH = c_right(i, j, k);
        const double LminH = u_right(i, j, k).x() - cH;
        const double LmaxH = u_right(i, j, k).x() + cH;

        const double f_Erad_T{erad_left(i - 1, j, k) * u_left(i - 1, j, k).x()};
        const double f_Erad_H{erad_right(i, j, k) * u_right(i, j, k).x()};

        // clang-format off
        erad_f(i, j, k) = numerical_algorithms::advect_conserved<double>(
          erad_left(i - 1, j, k), erad_right(i, j, k), f_Erad_T, f_Erad_H,
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
          one_over_dx_i[0] * (erad_f(i, j, k) - erad_f(i + 1, j, k));
      }; // forall
    }
    else if(face_axis == 1) {
      const auto mdpolicy_qcq = get_mdiota_policy(erad_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qcq) {

        auto [k, j, i] = kji;

        // min/max characteristic speeds on left
        const double cT = c_left(i, j - 1, k);
        const double LminT = u_left(i, j - 1, k).y() - cT;
        const double LmaxT = u_left(i, j - 1, k).y() + cT;

        // min/max characteristic speeds on right
        const double cH = c_right(i, j, k);
        const double LminH = u_right(i, j, k).y() - cH;
        const double LmaxH = u_right(i, j, k).y() + cH;

        const double f_Erad_T{erad_left(i, j - 1, k) * u_left(i, j - 1, k).y()};
        const double f_Erad_H{erad_right(i, j, k) * u_right(i, j, k).y()};

        // clang-format off
        erad_f(i, j, k) = numerical_algorithms::advect_conserved<double>(
          erad_left(i, j - 1, k), erad_right(i, j, k), f_Erad_T, f_Erad_H, LminT,
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
          one_over_dx_i[1] * (erad_f(i, j, k) - erad_f(i, j + 1, k));
      }; // forall
    }
    else {
      const auto mdpolicy_cqq = get_mdiota_policy(erad_f,
        m.template cells<ax::z, dm::corrector>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_cqq) {

        auto [k, j, i] = kji;

        // min/max characteristic speeds on left
        const double cT = c_left(i, j, k - 1);
        const double LminT = u_left(i, j, k - 1).z() - cT;
        const double LmaxT = u_left(i, j, k - 1).z() + cT;

        // min/max characteristic speeds on right
        const double cH = c_right(i, j, k);
        const double LminH = u_right(i, j, k).z() - cH;
        const double LmaxH = u_right(i, j, k).z() + cH;

        const double f_Erad_T{erad_left(i, j, k - 1) * u_left(i, j, k - 1).z()};
        const double f_Erad_H{erad_right(i, j, k) * u_right(i, j, k).z()};

        // clang-format off
        erad_f(i, j, k) = numerical_algorithms::advect_conserved<double>(
          erad_left(i, j, k - 1), erad_right(i, j, k), f_Erad_T, f_Erad_H,
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
          one_over_dx_i[2] * (erad_f(i, j, k) - erad_f(i, j, k + 1));
      }; // forall
    }
  }
}

} // namespace hard::tasks::rad

#endif
