#ifndef HARD_MODULE_HYDRO_INTERFACE_FLUXES_HH
#define HARD_MODULE_HYDRO_INTERFACE_FLUXES_HH

#include "../numerical_algorithms/riemann_solvers.hh"
#include <cstddef>

namespace hard::tasks::hydro {

template<std::size_t Dim>
void
compute_interface_fluxes(flecsi::exec::cpu s,
  std::size_t face_axis,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces<Dim>::accessor<ro, ro> r_face_a,
  typename faces_vec<Dim>::accessor<ro, ro> uFace_a,
  typename faces<Dim>::accessor<ro, ro> pFace_a,
  typename faces<Dim>::accessor<ro, ro> cFace_a,
  typename faces_vec<Dim>::accessor<ro, ro> ru_face_a,
  typename faces<Dim>::accessor<ro, ro> re_face_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> r_f_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ru_f_a,
  field<double>::accessor<wo, ro> re_f_a,
  // time derivative
  typename RK<Dim>::accessor<rw, na> rk_dt_a,
  typename single<vec<Dim>>::template accessor<ro> g_acc) noexcept {
  auto g = g_acc.get();
  auto [r_right, r_left] = faces<Dim>::mdcolex(m, r_face_a);
  auto [u_right, u_left] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [p_right, p_left] = faces<Dim>::mdcolex(m, pFace_a);
  auto [c_right, c_left] = faces<Dim>::mdcolex(m, cFace_a);
  auto [ru_right, ru_left] = faces_vec<Dim>::mdcolex(m, ru_face_a);
  auto [re_right, re_left] = faces<Dim>::mdcolex(m, re_face_a);

  auto r_f = m.template mdcolex<is::cells>(r_f_a);
  auto ru_f = m.template mdcolex<is::cells>(ru_f_a);
  auto re_f = m.template mdcolex<is::cells>(re_f_a);

  auto [dt_mass_density, dt_total_energy_density, dt_momentum_density] =
    RK<Dim>::mdcolex(m, rk_dt_a);

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
      const double pLeft_wave =
        p_left(i - 1) + 0.5 * r_left(i - 1) * g.x() * dx;
      const double pRight_wave = p_right(i) - 0.5 * r_right(i) * g.x() * dx;
      const double f_r_T{ru_left(i - 1).x()};
      const double f_r_H{ru_right(i).x()};
      const vec<1> f_ru_T{ru_left(i - 1).x() * u_left(i - 1).x() + pLeft_wave};
      const vec<1> f_ru_H{ru_right(i).x() * u_right(i).x() + pRight_wave};
      const double f_rE_T{(re_left(i - 1) + pLeft_wave) * u_left(i - 1).x()};
      const double f_rE_H{(re_right(i) + pRight_wave) * u_right(i).x()};

      // clang-format off
        r_f(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          r_left(i-1), r_left(i-1), u_left(i-1), re_left(i-1), pLeft_wave, c_left(i-1), f_r_T,
          r_right(i),   r_right(i),   u_right(i),   re_right(i),   pRight_wave, c_right(i), f_r_H,
          "rho");
        ru_f(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ru_left(i-1), r_left(i-1), u_left(i-1), re_left(i-1), pLeft_wave, c_left(i-1), f_ru_T,
          ru_right(i),   r_right(i),   u_right(i),   re_right(i),   pRight_wave, c_right(i), f_ru_H,
          "rhou");
        re_f(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          re_left(i-1), r_left(i-1), u_left(i-1), re_left(i-1), pLeft_wave, c_left(i-1), f_rE_T,
          re_right(i),   r_right(i),   u_right(i),   re_right(i),   pRight_wave, c_right(i), f_rE_H,
          "E");
      // clang-format on

    }; // forall

    // Store dF^x/dx into du_dt
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_mass_density(i) += one_over_dx_i[0] * (r_f(i) - r_f(i + 1));
      dt_momentum_density(i) += one_over_dx_i[0] * (ru_f(i) - ru_f(i + 1));
      dt_total_energy_density(i) += one_over_dx_i[0] * (re_f(i) - re_f(i + 1));
    }; // forall
  }
  else if constexpr(Dim == 2) {
#ifdef DEBUG
    if(face_axis > 1) {
      flog_fatal("2D executable cannot compute fluxes along Z axis");
    }
#endif
    if(face_axis == 0) {
      auto mdpolicy_qc = get_mdiota_policy(r_f,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(ji, mdpolicy_qc) {

        auto [j, i] = ji;

        // Fluxes from left and right state
        const double f_r_T{ru_left(i - 1, j).x()};
        const double f_r_H{ru_right(i, j).x()};
        const vec<2> f_ru_T{

          ru_left(i - 1, j).x() * u_left(i - 1, j).x() + p_left(i - 1, j),
          ru_left(i - 1, j).x() * u_left(i - 1, j).y()};
        const vec<2> f_ru_H{
          ru_right(i, j).x() * u_right(i, j).x() + p_right(i, j),
          ru_right(i, j).x() * u_right(i, j).y()};
        const double f_rE_T{
          (re_left(i - 1, j) + p_left(i - 1, j)) * u_left(i - 1, j).x()};
        const double f_rE_H{
          (re_right(i, j) + p_right(i, j)) * u_right(i, j).x()};

        // clang-format off
        const double pLeft_wave = p_left(i-1,j) - 0.5 * r_left(i-1,j) * g.x() / one_over_dx_i[0];
        const double pRight_wave = p_right(i,j) - 0.5 * r_right(i,j) * g.x() / one_over_dx_i[0];
        r_f(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          r_left(i-1,j), r_left(i-1,j), u_left(i-1,j), re_left(i-1,j), pLeft_wave, c_left(i-1,j), f_r_T,
          r_right(i,j),   r_right(i,j),   u_right(i,j),   re_right(i,j),   pRight_wave, c_right(i,j), f_r_H,
          "rho" );
        ru_f(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ru_left(i-1,j), r_left(i-1,j), u_left(i-1,j), re_left(i-1,j), pLeft_wave, c_left(i-1,j), f_ru_T,
          ru_right(i,j),   r_right(i,j),   u_right(i,j),   re_right(i,j),   pRight_wave, c_right(i,j), f_ru_H,
          "rhou" );
        re_f(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          re_left(i-1,j), r_left(i-1,j), u_left(i-1,j), re_left(i-1,j), pLeft_wave, c_left(i-1,j), f_rE_T,
          re_right(i,j),   r_right(i,j),   u_right(i,j),   re_right(i,j),   pRight_wave, c_right(i,j), f_rE_H,
          "E" );
        // clang-format on
      }; // forall

      // Store dF^x/dx into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;
        dt_mass_density(i, j) += one_over_dx_i[0] * (r_f(i, j) - r_f(i + 1, j));
        dt_momentum_density(i, j) +=
          one_over_dx_i[0] * (ru_f(i, j) - ru_f(i + 1, j));
        dt_total_energy_density(i, j) +=
          one_over_dx_i[0] * (re_f(i, j) - re_f(i + 1, j));

      }; // forall
    }
    else if(face_axis == 1) {
      auto mdpolicy_cq = get_mdiota_policy(r_f,
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_cq) {

        auto [j, i] = ji;

        // Fluxes from left and right state
        const double f_r_T{ru_left(i, j - 1).y()};
        const double f_r_H{ru_right(i, j).y()};
        const vec<2> f_ru_T{ru_left(i, j - 1).y() * u_left(i, j - 1).x(),
          ru_left(i, j - 1).y() * u_left(i, j - 1).y() + p_left(i, j - 1)};
        const vec<2> f_ru_H{ru_right(i, j).y() * u_right(i, j).x(),
          ru_right(i, j).y() * u_right(i, j).y() + p_right(i, j)};
        const double f_rE_T{
          (re_left(i, j - 1) + p_left(i, j - 1)) * u_left(i, j - 1).y()};
        const double f_rE_H{
          (re_right(i, j) + p_right(i, j)) * u_right(i, j).y()};

        // clang-format off
        const double pLeft_wave = p_left(i,j-1) - 0.5 * r_left(i,j-1) * g.y() / one_over_dx_i[1];
        const double pRight_wave = p_right(i,j) - 0.5 * r_right(i,j) * g.y() / one_over_dx_i[1];
        r_f(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          r_left(i,j-1), r_left(i,j-1), u_left(i,j-1), re_left(i,j-1), pLeft_wave, c_left(i,j-1), f_r_T,
          r_right(i,j),   r_right(i,j),   u_right(i,j),   re_right(i,j),   pRight_wave, c_right(i,j), f_r_H,
          "rho" );
        ru_f(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ru_left(i,j-1), r_left(i,j-1), u_left(i,j-1), re_left(i,j-1), pLeft_wave, c_left(i,j-1), f_ru_T,
          ru_right(i,j),   r_right(i,j),   u_right(i,j),   re_right(i,j),   pRight_wave, c_right(i,j), f_ru_H,
          "rhou" );
        re_f(i, j) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          re_left(i,j-1), r_left(i,j-1), u_left(i,j-1), re_left(i,j-1), pLeft_wave, c_left(i,j-1), f_rE_T,
          re_right(i,j),   r_right(i,j),   u_right(i,j),   re_right(i,j),   pRight_wave, c_right(i,j), f_rE_H,
          "E" );
        // clang-format on
      }; // forall

      // Store dF^y/dy into du_dt
      auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_qq) {
        auto [j, i] = ji;
        dt_mass_density(i, j) += one_over_dx_i[1] * (r_f(i, j) - r_f(i, j + 1));
        dt_momentum_density(i, j) +=
          one_over_dx_i[1] * (ru_f(i, j) - ru_f(i, j + 1));
        dt_total_energy_density(i, j) +=
          one_over_dx_i[1] * (re_f(i, j) - re_f(i, j + 1));

      }; // forall
    }
  }
  else { // Dim == 3

    if(face_axis == 0) {
      const auto mdpolicy_qqc = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(kji, mdpolicy_qqc) {

        auto [k, j, i] = kji;

        // Fluxes from left and right state
        const double f_r_T{ru_left(i - 1, j, k).x()};
        const double f_r_H{ru_right(i, j, k).x()};
        const vec<3> f_ru_T{ru_left(i - 1, j, k).x() * u_left(i - 1, j, k).x() +
                              p_left(i - 1, j, k),
          ru_left(i - 1, j, k).x() * u_left(i - 1, j, k).y(),
          ru_left(i - 1, j, k).x() * u_left(i - 1, j, k).z()};
        const vec<3> f_ru_H{
          ru_right(i, j, k).x() * u_right(i, j, k).x() + p_right(i, j, k),
          ru_right(i, j, k).x() * u_right(i, j, k).y(),
          ru_right(i, j, k).x() * u_right(i, j, k).z()};
        const double f_rE_T{(re_left(i - 1, j, k) + p_left(i - 1, j, k)) *
                            u_left(i - 1, j, k).x()};
        const double f_rE_H{
          (re_right(i, j, k) + p_right(i, j, k)) * u_right(i, j, k).x()};

        // Advect conserved quantities
        // clang-format off
        const double pLeft_wave = p_left(i-1,j,k) + 0.5 * r_left(i-1,j,k) * g.x() / one_over_dx_i[0];
        const double pRight_wave = p_right(i,j,k) - 0.5 * r_right(i,j,k) * g.x() / one_over_dx_i[0];
        r_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          r_left(i-1,j,k), r_left(i-1,j,k), u_left(i-1,j,k), re_left(i-1,j,k), pLeft_wave, c_left(i-1,j,k), f_r_T,
          r_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_r_H,
          "rho" );
        ru_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ru_left(i-1,j,k), r_left(i-1,j,k), u_left(i-1,j,k), re_left(i-1,j,k), pLeft_wave, c_left(i-1,j,k), f_ru_T,
          ru_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_ru_H,
          "rhou" );
        re_f(i, j,k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          re_left(i-1,j,k), r_left(i-1,j,k), u_left(i-1,j,k), re_left(i-1,j,k), pLeft_wave, c_left(i-1,j,k), f_rE_T,
          re_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_rE_H,
          "E" );
        // clang-format on

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
          one_over_dx_i[0] * (r_f(i, j, k) - r_f(i + 1, j, k));
        dt_momentum_density(i, j, k) +=
          one_over_dx_i[0] * (ru_f(i, j, k) - ru_f(i + 1, j, k));
        dt_total_energy_density(i, j, k) +=
          one_over_dx_i[0] * (re_f(i, j, k) - re_f(i + 1, j, k));
      }; // forall
    }
    else if(face_axis == 1) {
      const auto mdpolicy_qcq = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qcq) {

        auto [k, j, i] = kji;

        // Fluxes from left and right state
        const double f_r_T{ru_left(i, j - 1, k).y()};
        const double f_r_H{ru_right(i, j, k).y()};
        const vec<3> f_ru_T{ru_left(i, j - 1, k).y() * u_left(i, j - 1, k).x(),
          ru_left(i, j - 1, k).y() * u_left(i, j - 1, k).y() +
            p_left(i, j - 1, k),
          ru_left(i, j - 1, k).y() * u_left(i, j - 1, k).z()};
        const vec<3> f_ru_H{ru_right(i, j, k).y() * u_right(i, j, k).x(),
          ru_right(i, j, k).y() * u_right(i, j, k).y() + p_right(i, j, k),
          ru_right(i, j, k).y() * u_right(i, j, k).z()};
        const double f_rE_T{(re_left(i, j - 1, k) + p_left(i, j - 1, k)) *
                            u_left(i, j - 1, k).y()};
        const double f_rE_H{
          (re_right(i, j, k) + p_right(i, j, k)) * u_right(i, j, k).y()};

        // Advect conserved quantities

        // clang-format off
        const double pLeft_wave = p_left(i,j-1,k) + 0.5 * r_left(i,j-1,k) * g.y() / one_over_dx_i[1];
        const double pRight_wave = p_right(i,j,k) - 0.5 * r_right(i,j,k) * g.y() / one_over_dx_i[1];
        r_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          r_left(i,j-1,k), r_left(i,j-1,k), u_left(i,j-1,k), re_left(i,j-1,k), pLeft_wave, c_left(i,j-1,k), f_r_T,
          r_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_r_H,
          "rho" );
        ru_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ru_left(i,j-1,k), r_left(i,j-1,k), u_left(i,j-1,k), re_left(i,j-1,k), pLeft_wave, c_left(i,j-1,k), f_ru_T,
          ru_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_ru_H,
          "rhou" );
        re_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          re_left(i,j-1,k), r_left(i,j-1,k), u_left(i,j-1,k), re_left(i,j-1,k), pLeft_wave, c_left(i,j-1,k), f_rE_T,
          re_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_rE_H,
          "E" );
        // clang-format on
      }; // forall

      // Store dF^y/dy into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_mass_density(i, j, k) +=
          one_over_dx_i[1] * (r_f(i, j, k) - r_f(i, j + 1, k));
        dt_momentum_density(i, j, k) +=
          one_over_dx_i[1] * (ru_f(i, j, k) - ru_f(i, j + 1, k));
        dt_total_energy_density(i, j, k) +=
          one_over_dx_i[1] * (re_f(i, j, k) - re_f(i, j + 1, k));
      }; // forall
    }
    else {
      const auto mdpolicy_cqq = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::corrector>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_cqq) {

        auto [k, j, i] = kji;

        // Fluxes from left and right state
        const double f_r_T{ru_left(i, j, k - 1).z()};
        const double f_r_H{ru_right(i, j, k).z()};
        const vec<3> f_ru_T{ru_left(i, j, k - 1).z() * u_left(i, j, k - 1).x(),
          ru_left(i, j, k - 1).z() * u_left(i, j, k - 1).y(),
          ru_left(i, j, k - 1).z() * u_left(i, j, k - 1).z() +
            p_left(i, j, k - 1)};
        const vec<3> f_ru_H{ru_right(i, j, k).z() * u_right(i, j, k).x(),
          ru_right(i, j, k).z() * u_right(i, j, k).y(),
          ru_right(i, j, k).z() * u_right(i, j, k).z() + p_right(i, j, k)};
        const double f_rE_T{(re_left(i, j, k - 1) + p_left(i, j, k - 1)) *
                            u_left(i, j, k - 1).z()};
        const double f_rE_H{
          (re_right(i, j, k) + p_right(i, j, k)) * u_right(i, j, k).z()};

        // Advect conserved quantities

        // clang-format off
        const double pLeft_wave = p_left(i,j,k-1) + 0.5 * r_left(i,j,k-1) * g.z() / one_over_dx_i[2];
        const double pRight_wave = p_right(i,j,k) - 0.5 * r_right(i,j,k) * g.z() / one_over_dx_i[2];
        r_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          r_left(i,j,k-1), r_left(i,j,k-1), u_left(i,j,k-1), re_left(i,j,k-1), pLeft_wave, c_left(i,j,k-1), f_r_T,
          r_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_r_H,
          "rho" );
        ru_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis, 
          ru_left(i,j,k-1), r_left(i,j,k-1), u_left(i,j,k-1), re_left(i,j,k-1), pLeft_wave, c_left(i,j,k-1), f_ru_T,
          ru_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_ru_H,
          "rhou" );
        re_f(i, j, k) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis, 
          re_left(i,j,k-1), r_left(i,j,k-1), u_left(i,j,k-1), re_left(i,j,k-1), pLeft_wave, c_left(i,j,k-1), f_rE_T,
          re_right(i,j,k),   r_right(i,j,k),   u_right(i,j,k),   re_right(i,j,k),   pRight_wave, c_right(i,j,k), f_rE_H,
          "E" );
        // clang-format on
      }; // forall

      // Store dF^z/dz into du_dt
      const auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qqq) {
        auto [k, j, i] = kji;
        dt_mass_density(i, j, k) +=
          one_over_dx_i[2] * (r_f(i, j, k) - r_f(i, j, k + 1));
        dt_momentum_density(i, j, k) +=
          one_over_dx_i[2] * (ru_f(i, j, k) - ru_f(i, j, k + 1));
        dt_total_energy_density(i, j, k) +=
          one_over_dx_i[2] * (re_f(i, j, k) - re_f(i, j, k + 1));
      }; // forall
    }
  }
}

} // namespace hard::tasks::hydro

#endif // HARD_MODULE_HYDRO_INTERFACE_FLUXES_HH
