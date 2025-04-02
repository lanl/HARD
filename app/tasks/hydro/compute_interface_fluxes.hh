
#pragma once

#include "../../numerical_algorithms/riemann_solvers.hh"
#include "../../types.hh"
#include <cstddef>

namespace hard::tasks::hydro {

//
// Based on the reconstructed primitives and conservatives on cell interfaces,
// compute Riemann fluxes and store them into `*F` variables, then store the
// time derivative - (F^{i+1/2} - F^{i-1/2}) / dx into `dudt_fluxes`.
//
template<std::size_t Dim>
void
compute_interface_fluxes(std::size_t face_axis,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, ro> rTail_a,
  field<double>::accessor<wo, ro> rHead_a,
  typename field<vec<Dim>>::template accessor<wo, ro> uTail_a,
  typename field<vec<Dim>>::template accessor<wo, ro> uHead_a,
  field<double>::accessor<wo, ro> pTail_a,
  field<double>::accessor<wo, ro> pHead_a,
  field<double>::accessor<wo, ro> eTail_a,
  field<double>::accessor<wo, ro> eHead_a,
  field<double>::accessor<wo, ro> cTail_a,
  field<double>::accessor<wo, ro> cHead_a,
  field<double>::accessor<wo, ro>
#ifdef ENABLE_RADIATION
    EradTail_a
#endif
  ,
  field<double>::accessor<wo, ro>
#ifdef ENABLE_RADIATION
    EradHead_a
#endif
  ,
  typename field<vec<Dim>>::template accessor<wo, ro> ruTail_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ruHead_a,
  field<double>::accessor<wo, ro> rETail_a,
  field<double>::accessor<wo, ro> rEHead_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> rF_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ruF_a,
  field<double>::accessor<wo, ro> rEF_a,
  field<double>::accessor<wo, ro>
#ifdef ENABLE_RADIATION
    EradF_a
#endif
  ,
  // time derivative
  field<double>::accessor<rw, na> dt_mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> dt_momentum_density_a,
  field<double>::accessor<rw, na> dt_total_energy_density_a,
  field<double>::accessor<rw, na>
#ifdef ENABLE_RADIATION
    dt_radiation_energy_density_a
#endif
) {

  auto rTail = m.template mdcolex<is::cells>(rTail_a);
  auto rHead = m.template mdcolex<is::cells>(rHead_a);
  auto uTail = m.template mdcolex<is::cells>(uTail_a);
  auto uHead = m.template mdcolex<is::cells>(uHead_a);
  auto pTail = m.template mdcolex<is::cells>(pTail_a);
  auto pHead = m.template mdcolex<is::cells>(pHead_a);
  auto eTail = m.template mdcolex<is::cells>(eTail_a);
  auto eHead = m.template mdcolex<is::cells>(eHead_a);
  auto cTail = m.template mdcolex<is::cells>(cTail_a);
  auto cHead = m.template mdcolex<is::cells>(cHead_a);

#ifdef ENABLE_RADIATION
  auto EradTail = m.template mdcolex<is::cells>(EradTail_a);
  auto EradHead = m.template mdcolex<is::cells>(EradHead_a);
#endif
  auto ruTail = m.template mdcolex<is::cells>(ruTail_a);
  auto ruHead = m.template mdcolex<is::cells>(ruHead_a);
  auto rETail = m.template mdcolex<is::cells>(rETail_a);
  auto rEHead = m.template mdcolex<is::cells>(rEHead_a);
  auto rF = m.template mdcolex<is::cells>(rF_a);
  auto ruF = m.template mdcolex<is::cells>(ruF_a);
  auto rEF = m.template mdcolex<is::cells>(rEF_a);
#ifdef ENABLE_RADIATION
  auto EradF = m.template mdcolex<is::cells>(EradF_a);
#endif
  auto dt_mass_density = m.template mdcolex<is::cells>(dt_mass_density_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
#ifdef ENABLE_RADIATION
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);
#endif
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
    forall(i, (m.template cells<ax::x, dm::corrector>()), "compute_flux_1d") {
      // min/max characteristic speeds on left

      const double cT = cTail(i - 1);
      const double LminT = uTail(i - 1).x - cT;
      const double LmaxT = uTail(i - 1).x + cT;

      // min/max characteristic speeds on right
      const double cH = cHead(i);
      const double LminH = uHead(i).x - cH;
      const double LmaxH = uHead(i).x + cH;

      // Fluxes from left and right state
      const double f_r_T{ruTail(i - 1).x};
      const double f_r_H{ruHead(i).x};
      const vec<1> f_ru_T{ruTail(i - 1).x * uTail(i - 1).x + pTail(i - 1)};
      const vec<1> f_ru_H{ruHead(i).x * uHead(i).x + pHead(i)};
      const double f_rE_T{(rETail(i - 1) + pTail(i - 1)) * uTail(i - 1).x};
      const double f_rE_H{(rEHead(i) + pHead(i)) * uHead(i).x};

      // clang-format off
      const auto hll_fluxes = numerical_algorithms::hll_hydro(
        rTail(i - 1), rHead(i), ruTail(i - 1), ruHead(i), rETail(i - 1), rEHead(i),
        f_r_T, f_r_H, f_ru_T, f_ru_H, f_rE_T, f_rE_H,
        LminT, LmaxT, LminH, LmaxH);
      // clang-format on

      rF(i) = get<0>(hll_fluxes);
      ruF(i) = get<1>(hll_fluxes);
      rEF(i) = get<2>(hll_fluxes);

#ifdef ENABLE_RADIATION
      EradF(i) = numerical_algorithms::advect_Erad(
        EradTail(i - 1), EradHead(i), uTail(i - 1).x, uHead(i).x);
#endif
    }; // forall

    // Store dF^x/dx into du_dt
    forall(i, (m.template cells<ax::x, dm::quantities>()), "update_flux_1d") {
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

      forall(ji, mdpolicy_qc, "compute_flux_2dx") {

        auto [j, i] = ji;
        // min/max characteristic speeds on left
        const double cT = cTail(i - 1, j);
        const double LminT = uTail(i - 1, j).x - cT;
        const double LmaxT = uTail(i - 1, j).x + cT;

        // min/max characteristic speeds on right
        const double cH = cHead(i, j);
        const double LminH = uHead(i, j).x - cH;
        const double LmaxH = uHead(i, j).x + cH;

        // Fluxes from left and right state
        const double f_r_T{ruTail(i - 1, j).x};
        const double f_r_H{ruHead(i, j).x};
        const vec<2> f_ru_T{
          ruTail(i - 1, j).x * uTail(i - 1, j).x + pTail(i - 1, j),
          ruTail(i - 1, j).x * uTail(i - 1, j).y};
        const vec<2> f_ru_H{ruHead(i, j).x * uHead(i, j).x + pHead(i, j),
          ruHead(i, j).x * uHead(i, j).y};
        const double f_rE_T{
          (rETail(i - 1, j) + pTail(i - 1, j)) * uTail(i - 1, j).x};
        const double f_rE_H{(rEHead(i, j) + pHead(i, j)) * uHead(i, j).x};

        // clang-format off
          const auto hll_fluxes =
            numerical_algorithms::hll_hydro(
              rTail(i - 1, j), rHead(i, j), ruTail(i - 1, j), ruHead(i, j), rETail(i - 1, j), rEHead(i, j),
              f_r_T, f_r_H, f_ru_T, f_ru_H, f_rE_T, f_rE_H,
              LminT, LmaxT, LminH, LmaxH);
        // clang-format on

        rF(i, j) = get<0>(hll_fluxes);
        ruF(i, j) = get<1>(hll_fluxes);
        rEF(i, j) = get<2>(hll_fluxes);

#ifdef ENABLE_RADIATION
        EradF(i, j) = numerical_algorithms::advect_Erad(
          EradTail(i - 1, j), EradHead(i, j), uTail(i - 1, j).x, uHead(i, j).x);
#endif

      }; // forall

      // Store dF^x/dx into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      forall(ji, mdpolicy_qq, "update_flux_2dx") {
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

      forall(ji, mdpolicy_cq, "compute_flux_2dy") {

        auto [j, i] = ji;
        // min/max characteristic speeds on left
        const double cT = cTail(i, j - 1);
        const double LminT = uTail(i, j - 1).y - cT;
        const double LmaxT = uTail(i, j - 1).y + cT;

        // min/max characteristic speeds on right
        const double cH = cHead(i, j);
        const double LminH = uHead(i, j).y - cH;
        const double LmaxH = uHead(i, j).y + cH;

        // Fluxes from left and right state
        const double f_r_T{ruTail(i, j - 1).y};
        const double f_r_H{ruHead(i, j).y};
        const vec<2> f_ru_T{ruTail(i, j - 1).y * uTail(i, j - 1).x,
          ruTail(i, j - 1).y * uTail(i, j - 1).y + pTail(i, j - 1)};
        const vec<2> f_ru_H{ruHead(i, j).y * uHead(i, j).x,
          ruHead(i, j).y * uHead(i, j).y + pHead(i, j)};
        const double f_rE_T{
          (rETail(i, j - 1) + pTail(i, j - 1)) * uTail(i, j - 1).y};
        const double f_rE_H{(rEHead(i, j) + pHead(i, j)) * uHead(i, j).y};

        // clang-format off
          const auto hll_fluxes =
            numerical_algorithms::hll_hydro(
              rTail(i, j - 1), rHead(i, j), ruTail(i, j - 1), ruHead(i, j), rETail(i, j - 1), rEHead(i, j),
              f_r_T, f_r_H, f_ru_T, f_ru_H, f_rE_T, f_rE_H,
              LminT, LmaxT, LminH, LmaxH);
        // clang-format on

        rF(i, j) = get<0>(hll_fluxes);
        ruF(i, j) = get<1>(hll_fluxes);
        rEF(i, j) = get<2>(hll_fluxes);

#ifdef ENABLE_RADIATION
        EradF(i, j) = numerical_algorithms::advect_Erad(
          EradTail(i, j - 1), EradHead(i, j), uTail(i, j - 1).y, uHead(i, j).y);
#endif

      }; // forall

      // Store dF^y/dy into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      forall(ji, mdpolicy_qq, "update_flux_2dy") {
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

      forall(kji, mdpolicy_qqc, "compute_flux_3dx") {

        auto [k, j, i] = kji;
        // min/max characteristic speeds on left
        const double cT = cTail(i - 1, j, k);
        const double LminT = uTail(i - 1, j, k).x - cT;
        const double LmaxT = uTail(i - 1, j, k).x + cT;

        // min/max characteristic speeds on right
        const double cH = cHead(i, j, k);
        const double LminH = uHead(i, j, k).x - cH;
        const double LmaxH = uHead(i, j, k).x + cH;

        // Fluxes from left and right state
        const double f_r_T{ruTail(i - 1, j, k).x};
        const double f_r_H{ruHead(i, j, k).x};
        const vec<3> f_ru_T{
          ruTail(i - 1, j, k).x * uTail(i - 1, j, k).x + pTail(i - 1, j, k),
          ruTail(i - 1, j, k).x * uTail(i - 1, j, k).y,
          ruTail(i - 1, j, k).x * uTail(i - 1, j, k).z};
        const vec<3> f_ru_H{
          ruHead(i, j, k).x * uHead(i, j, k).x + pHead(i, j, k),
          ruHead(i, j, k).x * uHead(i, j, k).y,
          ruHead(i, j, k).x * uHead(i, j, k).z};
        const double f_rE_T{
          (rETail(i - 1, j, k) + pTail(i - 1, j, k)) * uTail(i - 1, j, k).x};
        const double f_rE_H{
          (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).x};

        // clang-format off
            const auto hll_fluxes =
              numerical_algorithms::hll_hydro(
                rTail(i - 1, j, k), rHead(i, j, k),
                ruTail(i - 1, j, k), ruHead(i, j, k),
                rETail(i - 1, j, k), rEHead(i, j, k),
                f_r_T, f_r_H, f_ru_T, f_ru_H, f_rE_T, f_rE_H,
                LminT, LmaxT, LminH, LmaxH);
        // clang-format on

        rF(i, j, k) = get<0>(hll_fluxes);
        ruF(i, j, k) = get<1>(hll_fluxes);
        rEF(i, j, k) = get<2>(hll_fluxes);

#ifdef ENABLE_RADIATION
        EradF(i, j, k) =
          numerical_algorithms::advect_Erad(EradTail(i - 1, j, k),
            EradHead(i, j, k),
            uTail(i - 1, j, k).x,
            uHead(i, j, k).x);
#endif

      }; // forall

      // Store dF^x/dx into du_dt
      // Store dF^x/dx into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      forall(kji, mdpolicy_qqq, "update_flux_3dx") {
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

      forall(kji, mdpolicy_qcq, "compute_flux_3dy") {

        auto [k, j, i] = kji;
        // min/max characteristic speeds on left
        const double cT = cTail(i, j - 1, k);
        const double LminT = uTail(i, j - 1, k).y - cT;
        const double LmaxT = uTail(i, j - 1, k).y + cT;

        // min/max characteristic speeds on right
        const double cH = cHead(i, j, k);
        const double LminH = uHead(i, j, k).y - cH;
        const double LmaxH = uHead(i, j, k).y + cH;

        // Fluxes from left and right state
        const double f_r_T{ruTail(i, j - 1, k).y};
        const double f_r_H{ruHead(i, j, k).y};
        const vec<3> f_ru_T{ruTail(i, j - 1, k).y * uTail(i, j - 1, k).x,
          ruTail(i, j - 1, k).y * uTail(i, j - 1, k).y + pTail(i, j - 1, k),
          ruTail(i, j - 1, k).y * uTail(i, j - 1, k).z};
        const vec<3> f_ru_H{ruHead(i, j, k).y * uHead(i, j, k).x,
          ruHead(i, j, k).y * uHead(i, j, k).y + pHead(i, j, k),
          ruHead(i, j, k).y * uHead(i, j, k).z};
        const double f_rE_T{
          (rETail(i, j - 1, k) + pTail(i, j - 1, k)) * uTail(i, j - 1, k).y};
        const double f_rE_H{
          (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).y};

        // clang-format off
            const auto hll_fluxes =
              numerical_algorithms::hll_hydro(
                rTail(i, j - 1, k), rHead(i, j, k),
                ruTail(i, j - 1, k), ruHead(i, j, k),
                rETail(i, j - 1, k), rEHead(i, j, k),
                f_r_T, f_r_H, f_ru_T, f_ru_H, f_rE_T, f_rE_H,
                LminT, LmaxT, LminH, LmaxH);
        // clang-format on

        rF(i, j, k) = get<0>(hll_fluxes);
        ruF(i, j, k) = get<1>(hll_fluxes);
        rEF(i, j, k) = get<2>(hll_fluxes);

#ifdef ENABLE_RADIATION
        EradF(i, j, k) =
          numerical_algorithms::advect_Erad(EradTail(i, j - 1, k),
            EradHead(i, j, k),
            uTail(i, j - 1, k).y,
            uHead(i, j, k).y);
#endif

      }; // forall

      // Store dF^y/dy into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      forall(kji, mdpolicy_qqq, "update_flux_3dy") {
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

      forall(kji, mdpolicy_cqq, "compute_flux_3dz") {

        auto [k, j, i] = kji;
        // min/max characteristic speeds on left
        const double cT = cTail(i, j, k - 1);
        const double LminT = uTail(i, j, k - 1).z - cT;
        const double LmaxT = uTail(i, j, k - 1).z + cT;

        // min/max characteristic speeds on right
        const double cH = cHead(i, j, k);
        const double LminH = uHead(i, j, k).z - cH;
        const double LmaxH = uHead(i, j, k).z + cH;

        // Fluxes from left and right state
        const double f_r_T{ruTail(i, j, k - 1).z};
        const double f_r_H{ruHead(i, j, k).z};
        const vec<3> f_ru_T{ruTail(i, j, k - 1).z * uTail(i, j, k - 1).x,
          ruTail(i, j, k - 1).z * uTail(i, j, k - 1).y,
          ruTail(i, j, k - 1).z * uTail(i, j, k - 1).z + pTail(i, j, k - 1)};
        const vec<3> f_ru_H{ruHead(i, j, k).z * uHead(i, j, k).x,
          ruHead(i, j, k).z * uHead(i, j, k).y,
          ruHead(i, j, k).z * uHead(i, j, k).z + pHead(i, j, k)};
        const double f_rE_T{
          (rETail(i, j, k - 1) + pTail(i, j, k - 1)) * uTail(i, j, k - 1).z};
        const double f_rE_H{
          (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).z};

        // clang-format off
            const auto hll_fluxes =
              numerical_algorithms::hll_hydro(
                rTail(i, j, k - 1), rHead(i, j, k),
                ruTail(i, j, k - 1), ruHead(i, j, k),
                rETail(i, j, k - 1), rEHead(i, j, k),
                f_r_T, f_r_H, f_ru_T, f_ru_H, f_rE_T, f_rE_H,
                LminT, LmaxT, LminH, LmaxH);
        // clang-format on

        rF(i, j, k) = get<0>(hll_fluxes);
        ruF(i, j, k) = get<1>(hll_fluxes);
        rEF(i, j, k) = get<2>(hll_fluxes);

#ifdef ENABLE_RADIATION
        EradF(i, j, k) =
          numerical_algorithms::advect_Erad(EradTail(i, j, k - 1),
            EradHead(i, j, k),
            uTail(i, j, k - 1).z,
            uHead(i, j, k).z);
#endif

      }; // forall

      // Store dF^z/dz into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      forall(kji, mdpolicy_qqq, "update_flux_3dz") {
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
