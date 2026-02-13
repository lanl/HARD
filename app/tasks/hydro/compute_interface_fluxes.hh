
#pragma once

#include "../../numerical_algorithms/riemann_solvers.hh"
#include "../../types.hh"
#include <cstddef>

namespace hard::tasks::hydro {

template<std::size_t Dim>
void
compute_interface_fluxes(flecsi::exec::cpu s,
  std::size_t face_axis,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces<Dim>::accessor<wo, ro> rFace_a,
  typename faces_vec<Dim>::accessor<wo, ro> uFace_a,
  typename faces<Dim>::accessor<wo, ro> pFace_a,
  typename faces<Dim>::accessor<wo, ro> cFace_a,
  typename faces<Dim>::accessor<wo, ro> EradFace_a,
  typename faces_vec<Dim>::accessor<wo, ro> ruFace_a,
  typename faces<Dim>::accessor<wo, ro> rEFace_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> rF_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ruF_a,
  field<double>::accessor<wo, ro> rEF_a,
  field<double>::accessor<wo, ro> EradF_a,
  // time derivative
  typename RK<Dim>::accessor<rw, na> rk_dt_a,
  typename single<vec<Dim>>::template accessor<ro> g_acc) noexcept {
  auto g = g_acc.get();
  auto [rRight, rLeft] = faces<Dim>::mdcolex(m, rFace_a);
  auto [uRight, uLeft] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [pRight, pLeft] = faces<Dim>::mdcolex(m, pFace_a);
  auto [cRight, cLeft] = faces<Dim>::mdcolex(m, cFace_a);
  auto [EradRight, EradLeft] = faces<Dim>::mdcolex(m, EradFace_a);
  auto [ruRight, ruLeft] = faces_vec<Dim>::mdcolex(m, ruFace_a);
  auto [rERight, rELeft] = faces<Dim>::mdcolex(m, rEFace_a);

  auto rF = m.template mdcolex<is::cells>(rF_a);
  auto ruF = m.template mdcolex<is::cells>(ruF_a);
  auto rEF = m.template mdcolex<is::cells>(rEF_a);
  auto EradF = m.template mdcolex<is::cells>(EradF_a);

  auto [dt_mass_density,
    dt_total_energy_density,
    dt_radiation_energy_density,
    dt_momentum_density] = RK<Dim>::mdcolex(m, rk_dt_a);

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
      const auto dx = m.template delta<ax::x>();

      // Fluxes from left and right state
      const double pLeft_wave = pLeft(i - 1) + 0.5 * rLeft(i - 1) * g.x() * dx;
      const double pRight_wave = pRight(i) - 0.5 * rRight(i) * g.x() * dx;
      const double f_r_T{ruLeft(i - 1).x()};
      const double f_r_H{ruRight(i).x()};
      const vec<1> f_ru_T{ruLeft(i - 1).x() * uLeft(i - 1).x() + pLeft_wave};
      const vec<1> f_ru_H{ruRight(i).x() * uRight(i).x() + pRight_wave};
      const double f_rE_T{(rELeft(i - 1) + pLeft_wave) * uLeft(i - 1).x()};
      const double f_rE_H{(rERight(i) + pRight_wave) * uRight(i).x()};

      // clang-format off
        rF(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rLeft(i-1), rLeft(i-1), uLeft(i-1), rELeft(i-1), pLeft_wave, cLeft(i-1), f_r_T,
          rRight(i),   rRight(i),   uRight(i),   rERight(i),   pRight_wave, cRight(i), f_r_H,
          "rho");
        ruF(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ruLeft(i-1), rLeft(i-1), uLeft(i-1), rELeft(i-1), pLeft_wave, cLeft(i-1), f_ru_T,
          ruRight(i),   rRight(i),   uRight(i),   rERight(i),   pRight_wave, cRight(i), f_ru_H,
          "rhou");
        rEF(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rELeft(i-1), rLeft(i-1), uLeft(i-1), rELeft(i-1), pLeft_wave, cLeft(i-1), f_rE_T,
          rERight(i),   rRight(i),   uRight(i),   rERight(i),   pRight_wave, cRight(i), f_rE_H,
          "E");
      // clang-format on

#ifdef ENABLE_RADIATION

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
#endif
    }; // forall

    // Store dF^x/dx into du_dt
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_mass_density(i) += one_over_dx_i[0] * (rF(i) - rF(i + 1));
      dt_momentum_density(i) += one_over_dx_i[0] * (ruF(i) - ruF(i + 1));
      dt_total_energy_density(i) += one_over_dx_i[0] * (rEF(i) - rEF(i + 1));
#ifdef ENABLE_RADIATION
      dt_radiation_energy_density(i) +=
        one_over_dx_i[0] * (EradF(i) - EradF(i + 1));
#endif
    }; // forall
  }
  else if constexpr(Dim == 2) {
#ifdef DEBUG
    if(face_axis > 1) {
      flog_fatal("2D executable cannot compute fluxes along Z axis");
    }
#endif
    if(face_axis == 0) {
      auto mdpolicy_qc = get_mdiota_policy(rF,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(ji, mdpolicy_qc) {

        auto [j, i] = ji;

        // Fluxes from left and right state
        const double f_r_T{ruLeft(i - 1, j).x()};
        const double f_r_H{ruRight(i, j).x()};
        const vec<2> f_ru_T{

          ruLeft(i - 1, j).x() * uLeft(i - 1, j).x() + pLeft(i - 1, j),
          ruLeft(i - 1, j).x() * uLeft(i - 1, j).y()};
        const vec<2> f_ru_H{ruRight(i, j).x() * uRight(i, j).x() + pRight(i, j),
          ruRight(i, j).x() * uRight(i, j).y()};
        const double f_rE_T{
          (rELeft(i - 1, j) + pLeft(i - 1, j)) * uLeft(i - 1, j).x()};
        const double f_rE_H{(rERight(i, j) + pRight(i, j)) * uRight(i, j).x()};

        // clang-format off
        const double pLeft_wave = pLeft(i-1,j) - 0.5 * rLeft(i-1,j) * g.x() / one_over_dx_i[0];
        const double pRight_wave = pRight(i,j) - 0.5 * rRight(i,j) * g.x() / one_over_dx_i[0];
        rF(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rLeft(i-1,j), rLeft(i-1,j), uLeft(i-1,j), rELeft(i-1,j), pLeft_wave, cLeft(i-1,j), f_r_T,
          rRight(i,j),   rRight(i,j),   uRight(i,j),   rERight(i,j),   pRight_wave, cRight(i,j), f_r_H,
          "rho" );
        ruF(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ruLeft(i-1,j), rLeft(i-1,j), uLeft(i-1,j), rELeft(i-1,j), pLeft_wave, cLeft(i-1,j), f_ru_T,
          ruRight(i,j),   rRight(i,j),   uRight(i,j),   rERight(i,j),   pRight_wave, cRight(i,j), f_ru_H,
          "rhou" );
        rEF(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rELeft(i-1,j), rLeft(i-1,j), uLeft(i-1,j), rELeft(i-1,j), pLeft_wave, cLeft(i-1,j), f_rE_T,
          rERight(i,j),   rRight(i,j),   uRight(i,j),   rERight(i,j),   pRight_wave, cRight(i,j), f_rE_H,
          "E" );
        // clang-format on

#ifdef ENABLE_RADIATION
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
#endif

      }; // forall

      // Store dF^x/dx into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;
        dt_mass_density(i, j) += one_over_dx_i[0] * (rF(i, j) - rF(i + 1, j));
        dt_momentum_density(i, j) +=
          one_over_dx_i[0] * (ruF(i, j) - ruF(i + 1, j));
        dt_total_energy_density(i, j) +=
          one_over_dx_i[0] * (rEF(i, j) - rEF(i + 1, j));
#ifdef ENABLE_RADIATION
        dt_radiation_energy_density(i, j) +=
          one_over_dx_i[0] * (EradF(i, j) - EradF(i + 1, j));
#endif

      }; // forall
    }
    else if(face_axis == 1) {
      auto mdpolicy_cq = get_mdiota_policy(rF,
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_cq) {

        auto [j, i] = ji;

        // Fluxes from left and right state
        const double f_r_T{ruLeft(i, j - 1).y()};
        const double f_r_H{ruRight(i, j).y()};
        const vec<2> f_ru_T{ruLeft(i, j - 1).y() * uLeft(i, j - 1).x(),
          ruLeft(i, j - 1).y() * uLeft(i, j - 1).y() + pLeft(i, j - 1)};
        const vec<2> f_ru_H{ruRight(i, j).y() * uRight(i, j).x(),
          ruRight(i, j).y() * uRight(i, j).y() + pRight(i, j)};
        const double f_rE_T{
          (rELeft(i, j - 1) + pLeft(i, j - 1)) * uLeft(i, j - 1).y()};
        const double f_rE_H{(rERight(i, j) + pRight(i, j)) * uRight(i, j).y()};

        // clang-format off
        const double pLeft_wave = pLeft(i,j-1) - 0.5 * rLeft(i,j-1) * g.y() / one_over_dx_i[1];
        const double pRight_wave = pRight(i,j) - 0.5 * rRight(i,j) * g.y() / one_over_dx_i[1];
        rF(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rLeft(i,j-1), rLeft(i,j-1), uLeft(i,j-1), rELeft(i,j-1), pLeft_wave, cLeft(i,j-1), f_r_T,
          rRight(i,j),   rRight(i,j),   uRight(i,j),   rERight(i,j),   pRight_wave, cRight(i,j), f_r_H,
          "rho" );
        ruF(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ruLeft(i,j-1), rLeft(i,j-1), uLeft(i,j-1), rELeft(i,j-1), pLeft_wave, cLeft(i,j-1), f_ru_T,
          ruRight(i,j),   rRight(i,j),   uRight(i,j),   rERight(i,j),   pRight_wave, cRight(i,j), f_ru_H,
          "rhou" );
        rEF(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rELeft(i,j-1), rLeft(i,j-1), uLeft(i,j-1), rELeft(i,j-1), pLeft_wave, cLeft(i,j-1), f_rE_T,
          rERight(i,j),   rRight(i,j),   uRight(i,j),   rERight(i,j),   pRight_wave, cRight(i,j), f_rE_H,
          "E" );
        // clang-format on

#ifdef ENABLE_RADIATION

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
#endif

      }; // forall

      // Store dF^y/dy into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;
        dt_mass_density(i, j) += one_over_dx_i[1] * (rF(i, j) - rF(i, j + 1));
        dt_momentum_density(i, j) +=
          one_over_dx_i[1] * (ruF(i, j) - ruF(i, j + 1));
        dt_total_energy_density(i, j) +=
          one_over_dx_i[1] * (rEF(i, j) - rEF(i, j + 1));
#ifdef ENABLE_RADIATION
        dt_radiation_energy_density(i, j) +=
          one_over_dx_i[1] * (EradF(i, j) - EradF(i, j + 1));
#endif

      }; // forall
    }
  }
  else { // Dim == 3

    if(face_axis == 0) {
      const auto mdpolicy_qqc = get_mdiota_policy(rF,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(kji, mdpolicy_qqc) {

        auto [k, j, i] = kji;

        // Fluxes from left and right state
        const double f_r_T{ruLeft(i - 1, j, k).x()};
        const double f_r_H{ruRight(i, j, k).x()};
        const vec<3> f_ru_T{
          ruLeft(i - 1, j, k).x() * uLeft(i - 1, j, k).x() + pLeft(i - 1, j, k),
          ruLeft(i - 1, j, k).x() * uLeft(i - 1, j, k).y(),
          ruLeft(i - 1, j, k).x() * uLeft(i - 1, j, k).z()};
        const vec<3> f_ru_H{
          ruRight(i, j, k).x() * uRight(i, j, k).x() + pRight(i, j, k),
          ruRight(i, j, k).x() * uRight(i, j, k).y(),
          ruRight(i, j, k).x() * uRight(i, j, k).z()};
        const double f_rE_T{
          (rELeft(i - 1, j, k) + pLeft(i - 1, j, k)) * uLeft(i - 1, j, k).x()};
        const double f_rE_H{
          (rERight(i, j, k) + pRight(i, j, k)) * uRight(i, j, k).x()};

        // Advect conserved quantities
        // clang-format off
        const double pLeft_wave = pLeft(i-1,j,k) + 0.5 * rLeft(i-1,j,k) * g.x() / one_over_dx_i[0];
        const double pRight_wave = pRight(i,j,k) - 0.5 * rRight(i,j,k) * g.x() / one_over_dx_i[0];
        rF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rLeft(i-1,j,k), rLeft(i-1,j,k), uLeft(i-1,j,k), rELeft(i-1,j,k), pLeft_wave, cLeft(i-1,j,k), f_r_T,
          rRight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_r_H,
          "rho" );
        ruF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ruLeft(i-1,j,k), rLeft(i-1,j,k), uLeft(i-1,j,k), rELeft(i-1,j,k), pLeft_wave, cLeft(i-1,j,k), f_ru_T,
          ruRight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_ru_H,
          "rhou" );
        rEF(i, j,k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rELeft(i-1,j,k), rLeft(i-1,j,k), uLeft(i-1,j,k), rELeft(i-1,j,k), pLeft_wave, cLeft(i-1,j,k), f_rE_T,
          rERight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_rE_H,
          "E" );
        // clang-format on

#ifdef ENABLE_RADIATION

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
#endif

      }; // forall

      // Store dF^x/dx into du_dt
      // Store dF^x/dx into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_mass_density(i, j, k) +=
          one_over_dx_i[0] * (rF(i, j, k) - rF(i + 1, j, k));
        dt_momentum_density(i, j, k) +=
          one_over_dx_i[0] * (ruF(i, j, k) - ruF(i + 1, j, k));
        dt_total_energy_density(i, j, k) +=
          one_over_dx_i[0] * (rEF(i, j, k) - rEF(i + 1, j, k));
#ifdef ENABLE_RADIATION
        dt_radiation_energy_density(i, j, k) +=
          one_over_dx_i[0] * (EradF(i, j, k) - EradF(i + 1, j, k));
#endif
      }; // forall
    }
    else if(face_axis == 1) {
      const auto mdpolicy_qcq = get_mdiota_policy(rF,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qcq) {

        auto [k, j, i] = kji;

        // Fluxes from left and right state
        const double f_r_T{ruLeft(i, j - 1, k).y()};
        const double f_r_H{ruRight(i, j, k).y()};
        const vec<3> f_ru_T{ruLeft(i, j - 1, k).y() * uLeft(i, j - 1, k).x(),
          ruLeft(i, j - 1, k).y() * uLeft(i, j - 1, k).y() + pLeft(i, j - 1, k),
          ruLeft(i, j - 1, k).y() * uLeft(i, j - 1, k).z()};
        const vec<3> f_ru_H{ruRight(i, j, k).y() * uRight(i, j, k).x(),
          ruRight(i, j, k).y() * uRight(i, j, k).y() + pRight(i, j, k),
          ruRight(i, j, k).y() * uRight(i, j, k).z()};
        const double f_rE_T{
          (rELeft(i, j - 1, k) + pLeft(i, j - 1, k)) * uLeft(i, j - 1, k).y()};
        const double f_rE_H{
          (rERight(i, j, k) + pRight(i, j, k)) * uRight(i, j, k).y()};

        // Advect conserved quantities

        // clang-format off
        const double pLeft_wave = pLeft(i,j-1,k) + 0.5 * rLeft(i,j-1,k) * g.y() / one_over_dx_i[1];
        const double pRight_wave = pRight(i,j,k) - 0.5 * rRight(i,j,k) * g.y() / one_over_dx_i[1];
        rF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rLeft(i,j-1,k), rLeft(i,j-1,k), uLeft(i,j-1,k), rELeft(i,j-1,k), pLeft_wave, cLeft(i,j-1,k), f_r_T,
          rRight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_r_H,
          "rho" );
        ruF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ruLeft(i,j-1,k), rLeft(i,j-1,k), uLeft(i,j-1,k), rELeft(i,j-1,k), pLeft_wave, cLeft(i,j-1,k), f_ru_T,
          ruRight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_ru_H,
          "rhou" );
        rEF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rELeft(i,j-1,k), rLeft(i,j-1,k), uLeft(i,j-1,k), rELeft(i,j-1,k), pLeft_wave, cLeft(i,j-1,k), f_rE_T,
          rERight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_rE_H,
          "E" );
        // clang-format on

#ifdef ENABLE_RADIATION

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
#endif

      }; // forall

      // Store dF^y/dy into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_mass_density(i, j, k) +=
          one_over_dx_i[1] * (rF(i, j, k) - rF(i, j + 1, k));
        dt_momentum_density(i, j, k) +=
          one_over_dx_i[1] * (ruF(i, j, k) - ruF(i, j + 1, k));
        dt_total_energy_density(i, j, k) +=
          one_over_dx_i[1] * (rEF(i, j, k) - rEF(i, j + 1, k));
#ifdef ENABLE_RADIATION
        dt_radiation_energy_density(i, j, k) +=
          one_over_dx_i[1] * (EradF(i, j, k) - EradF(i, j + 1, k));
#endif
      }; // forall
    }
    else {
      const auto mdpolicy_cqq = get_mdiota_policy(rF,
        m.template cells<ax::z, dm::corrector>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_cqq) {

        auto [k, j, i] = kji;

        // Fluxes from left and right state
        const double f_r_T{ruLeft(i, j, k - 1).z()};
        const double f_r_H{ruRight(i, j, k).z()};
        const vec<3> f_ru_T{ruLeft(i, j, k - 1).z() * uLeft(i, j, k - 1).x(),
          ruLeft(i, j, k - 1).z() * uLeft(i, j, k - 1).y(),
          ruLeft(i, j, k - 1).z() * uLeft(i, j, k - 1).z() +
            pLeft(i, j, k - 1)};
        const vec<3> f_ru_H{ruRight(i, j, k).z() * uRight(i, j, k).x(),
          ruRight(i, j, k).z() * uRight(i, j, k).y(),
          ruRight(i, j, k).z() * uRight(i, j, k).z() + pRight(i, j, k)};
        const double f_rE_T{
          (rELeft(i, j, k - 1) + pLeft(i, j, k - 1)) * uLeft(i, j, k - 1).z()};
        const double f_rE_H{
          (rERight(i, j, k) + pRight(i, j, k)) * uRight(i, j, k).z()};

        // Advect conserved quantities

        // clang-format off
        const double pLeft_wave = pLeft(i,j,k-1) + 0.5 * rLeft(i,j,k-1) * g.z() / one_over_dx_i[2];
        const double pRight_wave = pRight(i,j,k) - 0.5 * rRight(i,j,k) * g.z() / one_over_dx_i[2];
        rF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rLeft(i,j,k-1), rLeft(i,j,k-1), uLeft(i,j,k-1), rELeft(i,j,k-1), pLeft_wave, cLeft(i,j,k-1), f_r_T,
          rRight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_r_H,
          "rho" );
        ruF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ruLeft(i,j,k-1), rLeft(i,j,k-1), uLeft(i,j,k-1), rELeft(i,j,k-1), pLeft_wave, cLeft(i,j,k-1), f_ru_T,
          ruRight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_ru_H,
          "rhou" );
        rEF(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          rELeft(i,j,k-1), rLeft(i,j,k-1), uLeft(i,j,k-1), rELeft(i,j,k-1), pLeft_wave, cLeft(i,j,k-1), f_rE_T,
          rERight(i,j,k),   rRight(i,j,k),   uRight(i,j,k),   rERight(i,j,k),   pRight_wave, cRight(i,j,k), f_rE_H,
          "E" );
        // clang-format on

#ifdef ENABLE_RADIATION

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
#endif

      }; // forall

      // Store dF^z/dz into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_mass_density(i, j, k) +=
          one_over_dx_i[2] * (rF(i, j, k) - rF(i, j, k + 1));
        dt_momentum_density(i, j, k) +=
          one_over_dx_i[2] * (ruF(i, j, k) - ruF(i, j, k + 1));
        dt_total_energy_density(i, j, k) +=
          one_over_dx_i[2] * (rEF(i, j, k) - rEF(i, j, k + 1));
#ifdef ENABLE_RADIATION
        dt_radiation_energy_density(i, j, k) +=
          one_over_dx_i[2] * (EradF(i, j, k) - EradF(i, j, k + 1));
#endif
      }; // forall
    }
  }
}

} // namespace hard::tasks::hydro
