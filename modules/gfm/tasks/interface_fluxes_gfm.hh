#ifndef HARD_MODULE_GFM_INTERFACE_FLUXES_HH
#define HARD_MODULE_GFM_INTERFACE_FLUXES_HH

#include "../../hydro/numerical_algorithms/riemann_solvers.hh"
#include "../../rad/tasks/utils.hh"
#include "../modules/common/tasks/utils.hh"
#include "../modules/spec/eos.hh"
#include "../state.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace hard::tasks::gfm {

// Create artificial state for GFM cut-face.
// We need to provide complete and EOS-consistent state
// for the ghost side to HLLC (sound speed, total E den, mom).
template<std::size_t Dim>
void
compute_state_from_rho_u_p(const std::vector<eos::eos_wrapper> & eos_mm,
  int mat_id,
  double rho,
  const vec<Dim> & u,
  double p,
  double & e,
  double & c,
  vec<Dim> & ru,
  double & rE) noexcept {

  const double rho_f = std::max(rho, 1e-14);
  const double p_f = std::max(p, 1e-14);

  const auto mat = static_cast<std::size_t>(mat_id);
  const auto & eos = eos_mm[mat];

  // Getting internal energy
  e = hard::tasks::util::find_sie(eos, rho_f, p_f);
  // Compute sound speed
  c = eos.cRhoSie(rho_f, e);

  double uu = 0.0;
  for(std::size_t d = 0; d < Dim; ++d) {
    ru[d] = rho_f * u[d];
    uu += u[d] * u[d];
  }
  // Compute total energy density
  rE = rho_f * e + 0.5 * rho_f * uu;
}

template<std::size_t Dim>
void
compute_interface_fluxes_gfm(flecsi::exec::cpu s,
  std::size_t face_axis,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces<Dim>::template accessor<ro, ro> r_face_a,
  typename faces_vec<Dim>::template accessor<ro, ro> uFace_a,
  typename faces<Dim>::template accessor<ro, ro> pFace_a,
  typename faces<Dim>::template accessor<ro, ro> cFace_a,
  typename faces_vec<Dim>::template accessor<ro, ro> ru_face_a,
  typename faces<Dim>::template accessor<ro, ro> re_face_a,
  // GFM related
  typename field<int>::accessor<ro, ro> material_a,
  typename faces_vec<Dim>::template accessor<ro, ro> normal_f_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> r_f_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ru_f_a,
  field<double>::accessor<wo, ro> re_f_a,
  // time derivative
  typename RK<Dim>::template accessor<rw, na> rk_dt_a,
  typename single<vec<Dim>>::template accessor<ro> g_acc,
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {

  auto g = g_acc.get();

  auto [r_right, r_left] = faces<Dim>::mdcolex(m, r_face_a);
  auto [u_right, u_left] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [p_right, p_left] = faces<Dim>::mdcolex(m, pFace_a);
  auto [c_right, c_left] = faces<Dim>::mdcolex(m, cFace_a);
  auto [ru_right, ru_left] = faces_vec<Dim>::mdcolex(m, ru_face_a);
  auto [re_right, re_left] = faces<Dim>::mdcolex(m, re_face_a);
  auto [n_right, n_left] = faces_vec<Dim>::mdcolex(m, normal_f_a);

  auto mat_id = m.template mdcolex<is::cells>(material_a);

  auto r_f = m.template mdcolex<is::cells>(r_f_a);
  auto ru_f = m.template mdcolex<is::cells>(ru_f_a);
  auto re_f = m.template mdcolex<is::cells>(re_f_a);

  auto [dt_mass_density, dt_total_energy_density, dt_momentum_density] =
    RK<Dim>::mdcolex(m, rk_dt_a);

  using common::tasks::utils::get_mdiota_policy;

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

  // Compute axis normal to fill n_face if level set geometry
  // is not defined
  auto normalize_face_normal = [&](const vec<Dim> & nin) {
    double nn = 0.0;
    for(std::size_t d = 0; d < Dim; ++d)
      nn += nin[d] * nin[d];
    nn = std::sqrt(nn);

    vec<Dim> nout{};
    if(nn > 1e-14) {
      for(std::size_t d = 0; d < Dim; ++d)
        nout[d] = nin[d] / nn;
    }
    else {
      nout[face_axis] = 1.0;
    }
    return nout;
  };

  // ---------------------------------------------------------------------------
  // 1D
  // ---------------------------------------------------------------------------
  if constexpr(Dim == 1) {
    flog_assert(
      face_axis == 0, "1D executable can only compute fluxes along X axis");

    s.executor().forall(i, (m.template cells<ax::x, dm::corrector>())) {
      const auto dx = m.template delta<ax::x>();
      const auto inv_dx = one_over_dx_i[0];

      const int mL = mat_id(i - 1);
      const int mR = mat_id(i);

      const double pLeft_wave =
        p_left(i - 1) + 0.5 * r_left(i - 1) * g.x() * dx;
      const double pRight_wave = p_right(i) - 0.5 * r_right(i) * g.x() * dx;

      if(mL == mR) {
        const double f_r_T{ru_left(i - 1).x()};
        const double f_r_H{ru_right(i).x()};

        const vec<1> f_ru_T{
          ru_left(i - 1).x() * u_left(i - 1).x() + pLeft_wave};
        const vec<1> f_ru_H{ru_right(i).x() * u_right(i).x() + pRight_wave};

        const double f_rE_T{(re_left(i - 1) + pLeft_wave) * u_left(i - 1).x()};
        const double f_rE_H{(re_right(i) + pRight_wave) * u_right(i).x()};

        r_f(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
            r_left(i - 1),
            r_left(i - 1),
            u_left(i - 1),
            re_left(i - 1),
            pLeft_wave,
            c_left(i - 1),
            f_r_T,
            r_right(i),
            r_right(i),
            u_right(i),
            re_right(i),
            pRight_wave,
            c_right(i),
            f_r_H,
            numerical_algorithms::field_tag::rho);

        ru_f(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
            ru_left(i - 1),
            r_left(i - 1),
            u_left(i - 1),
            re_left(i - 1),
            pLeft_wave,
            c_left(i - 1),
            f_ru_T,
            ru_right(i),
            r_right(i),
            u_right(i),
            re_right(i),
            pRight_wave,
            c_right(i),
            f_ru_H,
            numerical_algorithms::field_tag::ru);

        re_f(i) =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
            re_left(i - 1),
            r_left(i - 1),
            u_left(i - 1),
            re_left(i - 1),
            pLeft_wave,
            c_left(i - 1),
            f_rE_T,
            re_right(i),
            r_right(i),
            u_right(i),
            re_right(i),
            pRight_wave,
            c_right(i),
            f_rE_H,
            numerical_algorithms::field_tag::re);

        dt_mass_density(i - 1) -= inv_dx * r_f(i);
        dt_momentum_density(i - 1) -= inv_dx * ru_f(i);
        dt_total_energy_density(i - 1) -= inv_dx * re_f(i);

        dt_mass_density(i) += inv_dx * r_f(i);
        dt_momentum_density(i) += inv_dx * ru_f(i);
        dt_total_energy_density(i) += inv_dx * re_f(i);
      }
      else {
        vec<1> n_if{};
        n_if[0] = 0.5 * (n_left(i - 1).x() + n_right(i).x());
        n_if = normalize_face_normal(n_if);

        // Real left/right states
        const double rhoL = std::max(r_left(i - 1), 1e-14);
        const vec<1> uL = u_left(i - 1);
        const double pL = std::max(pLeft_wave, 1e-14);
        const double cL = std::max(c_left(i - 1), 1e-14);
        const vec<1> ruL = ru_left(i - 1);
        const double reL = std::max(re_left(i - 1), 1e-14);

        const double rhoR = std::max(r_right(i), 1e-14);
        const vec<1> uR = u_right(i);
        const double pR = std::max(pRight_wave, 1e-14);
        const double cR = std::max(c_right(i), 1e-14);
        const vec<1> ruR = ru_right(i);
        const double reR = std::max(re_right(i), 1e-14);

        // Ghost for left owner: keep left rho, match right normal vel and
        // pressure
        const double unL = uL.x() * n_if.x();
        const double unR = uR.x() * n_if.x();

        vec<1> uRg_for_L{};
        uRg_for_L.x() = (uL.x() - unL * n_if.x()) + unR * n_if.x();

        const double rhoRg_L = rhoL;
        const double pRg_L = pR;

        double eRg_L{};
        double cRg_L{};
        vec<1> ruRg_L{};
        double reRg_L{};

        compute_state_from_rho_u_p<1>(
          eos_mm, mL, rhoRg_L, uRg_for_L, pRg_L, eRg_L, cRg_L, ruRg_L, reRg_L);

        const double f_r_T_L{ruL.x()};
        const double f_r_H_L{ruRg_L.x()};
        const vec<1> f_ru_T_L{ruL.x() * uL.x() + pL};
        const vec<1> f_ru_H_L{ruRg_L.x() * uRg_for_L.x() + pRg_L};
        const double f_rE_T_L{(reL + pL) * uL.x()};
        const double f_rE_H_L{(reRg_L + pRg_L) * uRg_for_L.x()};

        const double r_f_L =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
            rhoL,
            rhoL,
            uL,
            reL,
            pL,
            cL,
            f_r_T_L,
            rhoRg_L,
            rhoRg_L,
            uRg_for_L,
            reRg_L,
            pRg_L,
            cRg_L,
            f_r_H_L,
            numerical_algorithms::field_tag::rho);

        const vec<1> ru_f_L =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
            ruL,
            rhoL,
            uL,
            reL,
            pL,
            cL,
            f_ru_T_L,
            ruRg_L,
            rhoRg_L,
            uRg_for_L,
            reRg_L,
            pRg_L,
            cRg_L,
            f_ru_H_L,
            numerical_algorithms::field_tag::ru);

        const double re_f_L =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
            reL,
            rhoL,
            uL,
            reL,
            pL,
            cL,
            f_rE_T_L,
            reRg_L,
            rhoRg_L,
            uRg_for_L,
            reRg_L,
            pRg_L,
            cRg_L,
            f_rE_H_L,
            numerical_algorithms::field_tag::re);

        // Ghost for right owner
        vec<1> n_rev{};
        n_rev.x() = -n_if.x();

        const double unR_rev = uR.x() * n_rev.x();
        const double unL_rev = uL.x() * n_rev.x();

        vec<1> uLg_for_R{};
        uLg_for_R.x() = (uR.x() - unR_rev * n_rev.x()) + unL_rev * n_rev.x();

        const double rhoLg_R = rhoR;
        const double pLg_R = pL;
        double eLg_R{};
        double cLg_R{};
        vec<1> ruLg_R{};
        double reLg_R{};

        compute_state_from_rho_u_p<1>(
          eos_mm, mR, rhoLg_R, uLg_for_R, pLg_R, eLg_R, cLg_R, ruLg_R, reLg_R);

        const double f_r_T_R{ruLg_R.x()};
        const double f_r_H_R{ruR.x()};
        const vec<1> f_ru_T_R{ruLg_R.x() * uLg_for_R.x() + pLg_R};
        const vec<1> f_ru_H_R{ruR.x() * uR.x() + pR};
        const double f_rE_T_R{(reLg_R + pLg_R) * uLg_for_R.x()};
        const double f_rE_H_R{(reR + pR) * uR.x()};

        const double r_f_R =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
            rhoLg_R,
            rhoLg_R,
            uLg_for_R,
            reLg_R,
            pLg_R,
            cLg_R,
            f_r_T_R,
            rhoR,
            rhoR,
            uR,
            reR,
            pR,
            cR,
            f_r_H_R,
            numerical_algorithms::field_tag::rho);

        const vec<1> ru_f_R =
          numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
            ruLg_R,
            rhoLg_R,
            uLg_for_R,
            reLg_R,
            pLg_R,
            cLg_R,
            f_ru_T_R,
            ruR,
            rhoR,
            uR,
            reR,
            pR,
            cR,
            f_ru_H_R,
            numerical_algorithms::field_tag::ru);

        const double re_f_R =
          numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
            reLg_R,
            rhoLg_R,
            uLg_for_R,
            reLg_R,
            pLg_R,
            cLg_R,
            f_rE_T_R,
            reR,
            rhoR,
            uR,
            reR,
            pR,
            cR,
            f_rE_H_R,
            numerical_algorithms::field_tag::re);

        // store left-view flux
        r_f(i) = r_f_L;
        ru_f(i) = ru_f_L;
        re_f(i) = re_f_L;

        dt_mass_density(i - 1) -= inv_dx * r_f_L;
        dt_momentum_density(i - 1) -= inv_dx * ru_f_L;
        dt_total_energy_density(i - 1) -= inv_dx * re_f_L;

        dt_mass_density(i) += inv_dx * r_f_R;
        dt_momentum_density(i) += inv_dx * ru_f_R;
        dt_total_energy_density(i) += inv_dx * re_f_R;
      }
    };
  }

  // ---------------------------------------------------------------------------
  // 2D
  // ---------------------------------------------------------------------------
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
        const auto inv_dx = one_over_dx_i[0];

        const int mL = mat_id(i - 1, j);
        const int mR = mat_id(i, j);

        const double pLeft_wave =
          p_left(i - 1, j) + 0.5 * r_left(i - 1, j) * g.x() / one_over_dx_i[0];
        const double pRight_wave =
          p_right(i, j) - 0.5 * r_right(i, j) * g.x() / one_over_dx_i[0];

        // Same material face
        if(mL == mR) {
          const double f_r_T{ru_left(i - 1, j).x()};
          const double f_r_H{ru_right(i, j).x()};
          const vec<2> f_ru_T{
            ru_left(i - 1, j).x() * u_left(i - 1, j).x() + pLeft_wave,
            ru_left(i - 1, j).x() * u_left(i - 1, j).y()};
          const vec<2> f_ru_H{
            ru_right(i, j).x() * u_right(i, j).x() + pRight_wave,
            ru_right(i, j).x() * u_right(i, j).y()};
          const double f_rE_T{
            (re_left(i - 1, j) + pLeft_wave) * u_left(i - 1, j).x()};
          const double f_rE_H{
            (re_right(i, j) + pRight_wave) * u_right(i, j).x()};

          r_f(i, j) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              r_left(i - 1, j),
              r_left(i - 1, j),
              u_left(i - 1, j),
              re_left(i - 1, j),
              pLeft_wave,
              c_left(i - 1, j),
              f_r_T,
              r_right(i, j),
              r_right(i, j),
              u_right(i, j),
              re_right(i, j),
              pRight_wave,
              c_right(i, j),
              f_r_H,
              numerical_algorithms::field_tag::rho);

          ru_f(i, j) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ru_left(i - 1, j),
              r_left(i - 1, j),
              u_left(i - 1, j),
              re_left(i - 1, j),
              pLeft_wave,
              c_left(i - 1, j),
              f_ru_T,
              ru_right(i, j),
              r_right(i, j),
              u_right(i, j),
              re_right(i, j),
              pRight_wave,
              c_right(i, j),
              f_ru_H,
              numerical_algorithms::field_tag::ru);

          re_f(i, j) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              re_left(i - 1, j),
              r_left(i - 1, j),
              u_left(i - 1, j),
              re_left(i - 1, j),
              pLeft_wave,
              c_left(i - 1, j),
              f_rE_T,
              re_right(i, j),
              r_right(i, j),
              u_right(i, j),
              re_right(i, j),
              pRight_wave,
              c_right(i, j),
              f_rE_H,
              numerical_algorithms::field_tag::re);

          dt_mass_density(i - 1, j) -= inv_dx * r_f(i, j);
          dt_momentum_density(i - 1, j) -= inv_dx * ru_f(i, j);
          dt_total_energy_density(i - 1, j) -= inv_dx * re_f(i, j);

          dt_mass_density(i, j) += inv_dx * r_f(i, j);
          dt_momentum_density(i, j) += inv_dx * ru_f(i, j);
          dt_total_energy_density(i, j) += inv_dx * re_f(i, j);
        }
        // Cut face: GMF
        else {
          // Real left/right state
          const double rhoL = std::max(r_left(i - 1, j), 1e-14);
          const vec<2> uL = u_left(i - 1, j);
          const double pL = std::max(pLeft_wave, 1e-14);
          const double cL = std::max(c_left(i - 1, j), 1e-14);
          const vec<2> ruL = ru_left(i - 1, j);
          const double reL = std::max(re_left(i - 1, j), 1e-14);

          const double rhoR = std::max(r_right(i, j), 1e-14);
          const vec<2> uR = u_right(i, j);
          const double pR = std::max(pRight_wave, 1e-14);
          const double cR = std::max(c_right(i, j), 1e-14);
          const vec<2> ruR = ru_right(i, j);
          const double reR = std::max(re_right(i, j), 1e-14);

          // Interface normal from level set face field (average left/right
          // views)
          vec<2> n_if{};
          n_if.x() = 0.5 * (n_left(i - 1, j).x() + n_right(i, j).x());
          n_if.y() = 0.5 * (n_left(i - 1, j).y() + n_right(i, j).y());

          {
            const double nn =
              std::sqrt(n_if.x() * n_if.x() + n_if.y() * n_if.y());
            if(nn > 1e-14) {
              n_if.x() /= nn;
              n_if.y() /= nn;
            }
            else {
              n_if.x() = 1.0;
              n_if.y() = 0.0;
            }
          }

          // Left-owner ghost state
          // keep left density and tangential velocity
          // match right normal velocity and pressure
          const double unL = uL.x() * n_if.x() + uL.y() * n_if.y();
          const double unR = uR.x() * n_if.x() + uR.y() * n_if.y();

          vec<2> uRg_for_L{};
          uRg_for_L.x() = (uL.x() - unL * n_if.x()) + unR * n_if.x();
          uRg_for_L.y() = (uL.y() - unL * n_if.y()) + unR * n_if.y();

          const double rhoRg_L = rhoL;
          const double pRg_L = pR;
          double eRg_L{};
          double cRg_L{};
          vec<2> ruRg_L{};
          double reRg_L{};

          compute_state_from_rho_u_p<2>(eos_mm,
            mL,
            rhoRg_L,
            uRg_for_L,
            pRg_L,
            eRg_L,
            cRg_L,
            ruRg_L,
            reRg_L);

          const double f_r_T_L{ruL.x()};
          const double f_r_H_L{ruRg_L.x()};
          const vec<2> f_ru_T_L{ruL.x() * uL.x() + pL, ruL.x() * uL.y()};
          const vec<2> f_ru_H_L{
            ruRg_L.x() * uRg_for_L.x() + pRg_L, ruRg_L.x() * uRg_for_L.y()};
          const double f_rE_T_L{(reL + pL) * uL.x()};
          const double f_rE_H_L{(reRg_L + pRg_L) * uRg_for_L.x()};

          const double r_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_r_T_L,
              rhoRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_r_H_L,
              numerical_algorithms::field_tag::rho);

          const vec<2> ru_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_ru_T_L,
              ruRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_ru_H_L,
              numerical_algorithms::field_tag::ru);

          const double re_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_rE_T_L,
              reRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_rE_H_L,
              numerical_algorithms::field_tag::re);

          // Right owner ghost state
          // applying same idea, reversed normal
          vec<2> n_rev{};
          n_rev.x() = -n_if.x();
          n_rev.y() = -n_if.y();

          double unR_rev = uR.x() * n_rev.x() + uR.y() * n_rev.y();
          double unL_rev = uL.x() * n_rev.x() + uL.y() * n_rev.y();

          vec<2> uLg_for_R{};
          uLg_for_R.x() = (uR.x() - unR_rev * n_rev.x()) + unL_rev * n_rev.x();
          uLg_for_R.y() = (uR.y() - unR_rev * n_rev.y()) + unL_rev * n_rev.y();

          const double rhoLg_R = rhoR;
          const double pLg_R = pL;
          double eLg_R{};
          double cLg_R{};
          vec<2> ruLg_R{};
          double reLg_R{};

          compute_state_from_rho_u_p<2>(eos_mm,
            mR,
            rhoLg_R,
            uLg_for_R,
            pLg_R,
            eLg_R,
            cLg_R,
            ruLg_R,
            reLg_R);

          const double f_r_T_R{ruLg_R.x()};
          const double f_r_H_R{ruR.x()};
          const vec<2> f_ru_T_R{
            ruLg_R.x() * uLg_for_R.x() + pLg_R, ruLg_R.x() * uLg_for_R.y()};
          const vec<2> f_ru_H_R{ruR.x() * uR.x() + pR, ruR.x() * uR.y()};
          const double f_rE_T_R{(reLg_R + pLg_R) * uLg_for_R.x()};
          const double f_rE_H_R{(reR + pR) * uR.x()};

          const double r_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_r_T_R,
              rhoR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_r_H_R,
              numerical_algorithms::field_tag::rho);

          const vec<2> ru_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_ru_T_R,
              ruR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_ru_H_R,
              numerical_algorithms::field_tag::ru);

          const double re_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_rE_T_R,
              reR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_rE_H_R,
              numerical_algorithms::field_tag::re);

          r_f(i, j) = r_f_L;
          ru_f(i, j) = ru_f_L;
          re_f(i, j) = re_f_L;

          dt_mass_density(i - 1, j) -= inv_dx * r_f_L;
          dt_momentum_density(i - 1, j) -= inv_dx * ru_f_L;
          dt_total_energy_density(i - 1, j) -= inv_dx * re_f_L;

          dt_mass_density(i, j) += inv_dx * r_f_R;
          dt_momentum_density(i, j) += inv_dx * ru_f_R;
          dt_total_energy_density(i, j) += inv_dx * re_f_R;
        }
      };
    }
    else if(face_axis == 1) {
      auto mdpolicy_cq = get_mdiota_policy(r_f,
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(ji, mdpolicy_cq) {
        auto [j, i] = ji;
        const auto inv_dx = one_over_dx_i[1];

        const int mL = mat_id(i, j - 1);
        const int mR = mat_id(i, j);

        const double pLeft_wave =
          p_left(i, j - 1) + 0.5 * r_left(i, j - 1) * g.y() / one_over_dx_i[1];
        const double pRight_wave =
          p_right(i, j) - 0.5 * r_right(i, j) * g.y() / one_over_dx_i[1];
        // Same material face
        if(mL == mR) {
          const double f_r_T{ru_left(i, j - 1).y()};
          const double f_r_H{ru_right(i, j).y()};
          const vec<2> f_ru_T{ru_left(i, j - 1).y() * u_left(i, j - 1).x(),
            ru_left(i, j - 1).y() * u_left(i, j - 1).y() + pLeft_wave};
          const vec<2> f_ru_H{ru_right(i, j).y() * u_right(i, j).x(),
            ru_right(i, j).y() * u_right(i, j).y() + pRight_wave};
          const double f_rE_T{
            (re_left(i, j - 1) + pLeft_wave) * u_left(i, j - 1).y()};
          const double f_rE_H{
            (re_right(i, j) + pRight_wave) * u_right(i, j).y()};

          r_f(i, j) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              r_left(i, j - 1),
              r_left(i, j - 1),
              u_left(i, j - 1),
              re_left(i, j - 1),
              pLeft_wave,
              c_left(i, j - 1),
              f_r_T,
              r_right(i, j),
              r_right(i, j),
              u_right(i, j),
              re_right(i, j),
              pRight_wave,
              c_right(i, j),
              f_r_H,
              numerical_algorithms::field_tag::rho);

          ru_f(i, j) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ru_left(i, j - 1),
              r_left(i, j - 1),
              u_left(i, j - 1),
              re_left(i, j - 1),
              pLeft_wave,
              c_left(i, j - 1),
              f_ru_T,
              ru_right(i, j),
              r_right(i, j),
              u_right(i, j),
              re_right(i, j),
              pRight_wave,
              c_right(i, j),
              f_ru_H,
              numerical_algorithms::field_tag::ru);

          re_f(i, j) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              re_left(i, j - 1),
              r_left(i, j - 1),
              u_left(i, j - 1),
              re_left(i, j - 1),
              pLeft_wave,
              c_left(i, j - 1),
              f_rE_T,
              re_right(i, j),
              r_right(i, j),
              u_right(i, j),
              re_right(i, j),
              pRight_wave,
              c_right(i, j),
              f_rE_H,
              numerical_algorithms::field_tag::re);

          dt_mass_density(i, j - 1) -= inv_dx * r_f(i, j);
          dt_momentum_density(i, j - 1) -= inv_dx * ru_f(i, j);
          dt_total_energy_density(i, j - 1) -= inv_dx * re_f(i, j);

          dt_mass_density(i, j) += inv_dx * r_f(i, j);
          dt_momentum_density(i, j) += inv_dx * ru_f(i, j);
          dt_total_energy_density(i, j) += inv_dx * re_f(i, j);
        }
        // Cut face: GFM
        else {
          // Real left/right states
          const double rhoL = std::max(r_left(i, j - 1), 1e-14);
          const vec<2> uL = u_left(i, j - 1);
          const double pL = std::max(pLeft_wave, 1e-14);
          const double cL = std::max(c_left(i, j - 1), 1e-14);
          const vec<2> ruL = ru_left(i, j - 1);
          const double reL = std::max(re_left(i, j - 1), 1e-14);

          const double rhoR = std::max(r_right(i, j), 1e-14);
          const vec<2> uR = u_right(i, j);
          const double pR = std::max(pRight_wave, 1e-14);
          const double cR = std::max(c_right(i, j), 1e-14);
          const vec<2> ruR = ru_right(i, j);
          const double reR = std::max(re_right(i, j), 1e-14);

          // Interface normal from level set face field (average left/right
          // views)
          vec<2> n_if{};
          n_if.x() = 0.5 * (n_left(i, j - 1).x() + n_right(i, j).x());
          n_if.y() = 0.5 * (n_left(i, j - 1).y() + n_right(i, j).y());
          {
            const double nn =
              std::sqrt(n_if.x() * n_if.x() + n_if.y() * n_if.y());
            if(nn > 1e-14) {
              n_if.x() /= nn;
              n_if.y() /= nn;
            }
            else {
              n_if.x() = 0.0;
              n_if.y() = 1.0;
            }
          }
          // Left-owner ghost state
          const double unL = uL.x() * n_if.x() + uL.y() * n_if.y();
          const double unR = uR.x() * n_if.x() + uR.y() * n_if.y();

          vec<2> uRg_for_L{};
          uRg_for_L.x() = (uL.x() - unL * n_if.x()) + unR * n_if.x();
          uRg_for_L.y() = (uL.y() - unL * n_if.y()) + unR * n_if.y();

          const double rhoRg_L = rhoL;
          const double pRg_L = pR;
          double eRg_L{};
          double cRg_L{};
          vec<2> ruRg_L{};
          double reRg_L{};

          compute_state_from_rho_u_p<2>(eos_mm,
            mL,
            rhoRg_L,
            uRg_for_L,
            pRg_L,
            eRg_L,
            cRg_L,
            ruRg_L,
            reRg_L);

          const double f_r_T_L{ruL.y()};
          const double f_r_H_L{ruRg_L.y()};
          const vec<2> f_ru_T_L{ruL.y() * uL.x(), ruL.y() * uL.y() + pL};
          const vec<2> f_ru_H_L{
            ruRg_L.y() * uRg_for_L.x(), ruRg_L.y() * uRg_for_L.y() + pRg_L};
          const double f_rE_T_L{(reL + pL) * uL.y()};
          const double f_rE_H_L{(reRg_L + pRg_L) * uRg_for_L.y()};

          const double r_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_r_T_L,
              rhoRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_r_H_L,
              numerical_algorithms::field_tag::rho);

          const vec<2> ru_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_ru_T_L,
              ruRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_ru_H_L,
              numerical_algorithms::field_tag::ru);

          const double re_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_rE_T_L,
              reRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_rE_H_L,
              numerical_algorithms::field_tag::re);

          vec<2> n_rev{};
          n_rev.x() = -n_if.x();
          n_rev.y() = -n_if.y();

          double unR_rev = uR.x() * n_rev.x() + uR.y() * n_rev.y();
          double unL_rev = uL.x() * n_rev.x() + uL.y() * n_rev.y();

          vec<2> uLg_for_R{};
          uLg_for_R.x() = (uR.x() - unR_rev * n_rev.x()) + unL_rev * n_rev.x();
          uLg_for_R.y() = (uR.y() - unR_rev * n_rev.y()) + unL_rev * n_rev.y();

          const double rhoLg_R = rhoR;
          const double pLg_R = pL;
          double eLg_R{};
          double cLg_R{};
          vec<2> ruLg_R{};
          double reLg_R{};

          compute_state_from_rho_u_p<2>(eos_mm,
            mR,
            rhoLg_R,
            uLg_for_R,
            pLg_R,
            eLg_R,
            cLg_R,
            ruLg_R,
            reLg_R);

          const double f_r_T_R{ruLg_R.y()};
          const double f_r_H_R{ruR.y()};
          const vec<2> f_ru_T_R{
            ruLg_R.y() * uLg_for_R.x(), ruLg_R.y() * uLg_for_R.y() + pLg_R};
          const vec<2> f_ru_H_R{ruR.y() * uR.x(), ruR.y() * uR.y() + pR};
          const double f_rE_T_R{(reLg_R + pLg_R) * uLg_for_R.y()};
          const double f_rE_H_R{(reR + pR) * uR.y()};

          const double r_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_r_T_R,
              rhoR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_r_H_R,
              numerical_algorithms::field_tag::rho);

          const vec<2> ru_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_ru_T_R,
              ruR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_ru_H_R,
              numerical_algorithms::field_tag::ru);

          const double re_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_rE_T_R,
              reR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_rE_H_R,
              numerical_algorithms::field_tag::re);

          r_f(i, j) = r_f_L;
          ru_f(i, j) = ru_f_L;
          re_f(i, j) = re_f_L;

          dt_mass_density(i, j - 1) -= inv_dx * r_f_L;
          dt_momentum_density(i, j - 1) -= inv_dx * ru_f_L;
          dt_total_energy_density(i, j - 1) -= inv_dx * re_f_L;

          dt_mass_density(i, j) += inv_dx * r_f_R;
          dt_momentum_density(i, j) += inv_dx * ru_f_R;
          dt_total_energy_density(i, j) += inv_dx * re_f_R;
        }
      };
    }
  }

  // ---------------------------------------------------------------------------
  // 3D
  // ---------------------------------------------------------------------------
  else {
    if(face_axis == 0) {
      const auto mdpolicy_qqc = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      s.executor().forall(kji, mdpolicy_qqc) {
        auto [k, j, i] = kji;
        const auto inv_dx = one_over_dx_i[0];

        const int mL = mat_id(i - 1, j, k);
        const int mR = mat_id(i, j, k);

        const double pLeft_wave =
          p_left(i - 1, j, k) +
          0.5 * r_left(i - 1, j, k) * g.x() / one_over_dx_i[0];
        const double pRight_wave =
          p_right(i, j, k) - 0.5 * r_right(i, j, k) * g.x() / one_over_dx_i[0];

        if(mL == mR) {
          const double f_r_T{ru_left(i - 1, j, k).x()};
          const double f_r_H{ru_right(i, j, k).x()};
          const vec<3> f_ru_T{
            ru_left(i - 1, j, k).x() * u_left(i - 1, j, k).x() + pLeft_wave,
            ru_left(i - 1, j, k).x() * u_left(i - 1, j, k).y(),
            ru_left(i - 1, j, k).x() * u_left(i - 1, j, k).z()};
          const vec<3> f_ru_H{
            ru_right(i, j, k).x() * u_right(i, j, k).x() + pRight_wave,
            ru_right(i, j, k).x() * u_right(i, j, k).y(),
            ru_right(i, j, k).x() * u_right(i, j, k).z()};
          const double f_rE_T{
            (re_left(i - 1, j, k) + pLeft_wave) * u_left(i - 1, j, k).x()};
          const double f_rE_H{
            (re_right(i, j, k) + pRight_wave) * u_right(i, j, k).x()};

          r_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              r_left(i - 1, j, k),
              r_left(i - 1, j, k),
              u_left(i - 1, j, k),
              re_left(i - 1, j, k),
              pLeft_wave,
              c_left(i - 1, j, k),
              f_r_T,
              r_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_r_H,
              numerical_algorithms::field_tag::rho);

          ru_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ru_left(i - 1, j, k),
              r_left(i - 1, j, k),
              u_left(i - 1, j, k),
              re_left(i - 1, j, k),
              pLeft_wave,
              c_left(i - 1, j, k),
              f_ru_T,
              ru_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_ru_H,
              numerical_algorithms::field_tag::ru);

          re_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              re_left(i - 1, j, k),
              r_left(i - 1, j, k),
              u_left(i - 1, j, k),
              re_left(i - 1, j, k),
              pLeft_wave,
              c_left(i - 1, j, k),
              f_rE_T,
              re_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_rE_H,
              numerical_algorithms::field_tag::re);

          dt_mass_density(i - 1, j, k) -= inv_dx * r_f(i, j, k);
          dt_momentum_density(i - 1, j, k) -= inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i - 1, j, k) -= inv_dx * re_f(i, j, k);

          dt_mass_density(i, j, k) += inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k) += inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k) += inv_dx * re_f(i, j, k);
        }
        else {
          const double rhoL = std::max(r_left(i - 1, j, k), 1e-14);
          const vec<3> uL = u_left(i - 1, j, k);
          const double pL = std::max(pLeft_wave, 1e-14);
          const double cL = std::max(c_left(i - 1, j, k), 1e-14);
          const vec<3> ruL = ru_left(i - 1, j, k);
          const double reL = std::max(re_left(i - 1, j, k), 1e-14);

          const double rhoR = std::max(r_right(i, j, k), 1e-14);
          const vec<3> uR = u_right(i, j, k);
          const double pR = std::max(pRight_wave, 1e-14);
          const double cR = std::max(c_right(i, j, k), 1e-14);
          const vec<3> ruR = ru_right(i, j, k);
          const double reR = std::max(re_right(i, j, k), 1e-14);

          vec<3> n_if{};
          n_if.x() = 0.5 * (n_left(i - 1, j, k).x() + n_right(i, j, k).x());
          n_if.y() = 0.5 * (n_left(i - 1, j, k).y() + n_right(i, j, k).y());
          n_if.z() = 0.5 * (n_left(i - 1, j, k).z() + n_right(i, j, k).z());
          {
            const double nn = std::sqrt(
              n_if.x() * n_if.x() + n_if.y() * n_if.y() + n_if.z() * n_if.z());
            if(nn > 1e-14) {
              n_if.x() /= nn;
              n_if.y() /= nn;
              n_if.z() /= nn;
            }
            else {
              n_if.x() = 1.0;
              n_if.y() = 0.0;
              n_if.y() = 0.0;
            }
          }
          double unL =
            uL.x() * n_if.x() + uL.y() * n_if.y() + uL.z() * n_if.z();
          double unR =
            uR.x() * n_if.x() + uR.y() * n_if.y() + uR.z() * n_if.z();

          vec<3> uRg_for_L{};
          uRg_for_L.x() = (uL.x() - unL * n_if.x()) + unR * n_if.x();
          uRg_for_L.y() = (uL.y() - unL * n_if.y()) + unR * n_if.y();
          uRg_for_L.z() = (uL.z() - unL * n_if.z()) + unR * n_if.z();

          const double rhoRg_L = rhoL;
          const double pRg_L = pR;
          double eRg_L{};
          double cRg_L{};
          vec<3> ruRg_L{};
          double reRg_L{};

          compute_state_from_rho_u_p<3>(eos_mm,
            mL,
            rhoRg_L,
            uRg_for_L,
            pRg_L,
            eRg_L,
            cRg_L,
            ruRg_L,
            reRg_L);

          const double f_r_T_L{ruL.x()};
          const double f_r_H_L{ruRg_L.x()};
          const vec<3> f_ru_T_L{
            ruL.x() * uL.x() + pL, ruL.x() * uL.y(), ruL.x() * uL.z()};
          const vec<3> f_ru_H_L{ruRg_L.x() * uRg_for_L.x() + pRg_L,
            ruRg_L.x() * uRg_for_L.y(),
            ruRg_L.x() * uRg_for_L.z()};
          const double f_rE_T_L{(reL + pL) * uL.x()};
          const double f_rE_H_L{(reRg_L + pRg_L) * uRg_for_L.x()};

          const double r_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_r_T_L,
              rhoRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_r_H_L,
              numerical_algorithms::field_tag::rho);

          const vec<3> ru_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_ru_T_L,
              ruRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_ru_H_L,
              numerical_algorithms::field_tag::ru);

          const double re_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_rE_T_L,
              reRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_rE_H_L,
              numerical_algorithms::field_tag::re);

          vec<3> n_rev{};
          n_rev.x() = -n_if.x();
          n_rev.y() = -n_if.y();
          n_rev.z() = -n_if.z();

          double unR_rev =
            uR.x() * n_rev.x() + uR.y() * n_rev.y() + uR.z() * n_rev.z();
          double unL_rev =
            uL.x() * n_rev.x() + uL.y() * n_rev.y() + uL.z() * n_rev.z();

          vec<3> uLg_for_R{};
          uLg_for_R.x() = (uR.x() - unR_rev * n_rev.x()) + unL_rev * n_rev.x();
          uLg_for_R.y() = (uR.y() - unR_rev * n_rev.y()) + unL_rev * n_rev.y();
          uLg_for_R.z() = (uR.z() - unR_rev * n_rev.z()) + unL_rev * n_rev.z();

          const double rhoLg_R = rhoR;
          const double pLg_R = pL;
          double eLg_R{};
          double cLg_R{};
          vec<3> ruLg_R{};
          double reLg_R{};

          compute_state_from_rho_u_p<3>(eos_mm,
            mR,
            rhoLg_R,
            uLg_for_R,
            pLg_R,
            eLg_R,
            cLg_R,
            ruLg_R,
            reLg_R);

          const double f_r_T_R{ruLg_R.x()};
          const double f_r_H_R{ruR.x()};
          const vec<3> f_ru_T_R{ruLg_R.x() * uLg_for_R.x() + pLg_R,
            ruLg_R.x() * uLg_for_R.y(),
            ruLg_R.x() * uLg_for_R.z()};
          const vec<3> f_ru_H_R{
            ruR.x() * uR.x() + pR, ruR.x() * uR.y(), ruR.x() * uR.z()};
          const double f_rE_T_R{(reLg_R + pLg_R) * uLg_for_R.x()};
          const double f_rE_H_R{(reR + pR) * uR.x()};

          const double r_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_r_T_R,
              rhoR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_r_H_R,
              numerical_algorithms::field_tag::rho);

          const vec<3> ru_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_ru_T_R,
              ruR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_ru_H_R,
              numerical_algorithms::field_tag::ru);

          const double re_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_rE_T_R,
              reR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_rE_H_R,
              numerical_algorithms::field_tag::re);

          r_f(i, j, k) = r_f_L;
          ru_f(i, j, k) = ru_f_L;
          re_f(i, j, k) = re_f_L;

          dt_mass_density(i - 1, j, k) -= inv_dx * r_f_L;
          dt_momentum_density(i - 1, j, k) -= inv_dx * ru_f_L;
          dt_total_energy_density(i - 1, j, k) -= inv_dx * re_f_L;

          dt_mass_density(i, j, k) += inv_dx * r_f_R;
          dt_momentum_density(i, j, k) += inv_dx * ru_f_R;
          dt_total_energy_density(i, j, k) += inv_dx * re_f_R;
        }
      };
    }
    else if(face_axis == 1) {
      const auto mdpolicy_qcq = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_qcq) {
        auto [k, j, i] = kji;
        const auto inv_dx = one_over_dx_i[1];

        const int mL = mat_id(i, j - 1, k);
        const int mR = mat_id(i, j, k);

        const double pLeft_wave =
          p_left(i, j - 1, k) +
          0.5 * r_left(i, j - 1, k) * g.y() / one_over_dx_i[1];
        const double pRight_wave =
          p_right(i, j, k) - 0.5 * r_right(i, j, k) * g.y() / one_over_dx_i[1];

        if(mL == mR) {
          const double f_r_T{ru_left(i, j - 1, k).y()};
          const double f_r_H{ru_right(i, j, k).y()};
          const vec<3> f_ru_T{
            ru_left(i, j - 1, k).y() * u_left(i, j - 1, k).x(),
            ru_left(i, j - 1, k).y() * u_left(i, j - 1, k).y() + pLeft_wave,
            ru_left(i, j - 1, k).y() * u_left(i, j - 1, k).z()};
          const vec<3> f_ru_H{ru_right(i, j, k).y() * u_right(i, j, k).x(),
            ru_right(i, j, k).y() * u_right(i, j, k).y() + pRight_wave,
            ru_right(i, j, k).y() * u_right(i, j, k).z()};
          const double f_rE_T{
            (re_left(i, j - 1, k) + pLeft_wave) * u_left(i, j - 1, k).y()};
          const double f_rE_H{
            (re_right(i, j, k) + pRight_wave) * u_right(i, j, k).y()};

          r_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              r_left(i, j - 1, k),
              r_left(i, j - 1, k),
              u_left(i, j - 1, k),
              re_left(i, j - 1, k),
              pLeft_wave,
              c_left(i, j - 1, k),
              f_r_T,
              r_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_r_H,
              numerical_algorithms::field_tag::rho);

          ru_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ru_left(i, j - 1, k),
              r_left(i, j - 1, k),
              u_left(i, j - 1, k),
              re_left(i, j - 1, k),
              pLeft_wave,
              c_left(i, j - 1, k),
              f_ru_T,
              ru_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_ru_H,
              numerical_algorithms::field_tag::ru);

          re_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              re_left(i, j - 1, k),
              r_left(i, j - 1, k),
              u_left(i, j - 1, k),
              re_left(i, j - 1, k),
              pLeft_wave,
              c_left(i, j - 1, k),
              f_rE_T,
              re_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_rE_H,
              numerical_algorithms::field_tag::re);

          dt_mass_density(i, j - 1, k) -= inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j - 1, k) -= inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j - 1, k) -= inv_dx * re_f(i, j, k);

          dt_mass_density(i, j, k) += inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k) += inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k) += inv_dx * re_f(i, j, k);
        }
        else {
          const double rhoL = std::max(r_left(i, j - 1, k), 1e-14);
          const vec<3> uL = u_left(i, j - 1, k);
          const double pL = std::max(pLeft_wave, 1e-14);
          const double cL = std::max(c_left(i, j - 1, k), 1e-14);
          const vec<3> ruL = ru_left(i, j - 1, k);
          const double reL = std::max(re_left(i, j - 1, k), 1e-14);

          const double rhoR = std::max(r_right(i, j, k), 1e-14);
          const vec<3> uR = u_right(i, j, k);
          const double pR = std::max(pRight_wave, 1e-14);
          const double cR = std::max(c_right(i, j, k), 1e-14);
          const vec<3> ruR = ru_right(i, j, k);
          const double reR = std::max(re_right(i, j, k), 1e-14);

          vec<3> n_if{};
          n_if.x() = 0.5 * (n_left(i, j - 1, k).x() + n_right(i, j, k).x());
          n_if.y() = 0.5 * (n_left(i, j - 1, k).y() + n_right(i, j, k).y());
          n_if.z() = 0.5 * (n_left(i, j - 1, k).z() + n_right(i, j, k).z());
          {
            const double nn = std::sqrt(
              n_if.x() * n_if.x() + n_if.y() * n_if.y() + n_if.z() * n_if.z());
            if(nn > 1e-14) {
              n_if.x() /= nn;
              n_if.y() /= nn;
              n_if.z() /= nn;
            }
            else {
              n_if.x() = 0.0;
              n_if.y() = 1.0;
              n_if.y() = 0.0;
            }
          }
          double unL =
            uL.x() * n_if.x() + uL.y() * n_if.y() + uL.z() * n_if.z();
          double unR =
            uR.x() * n_if.x() + uR.y() * n_if.y() + uR.z() * n_if.z();

          vec<3> uRg_for_L{};
          uRg_for_L.x() = (uL.x() - unL * n_if.x()) + unR * n_if.x();
          uRg_for_L.y() = (uL.y() - unL * n_if.y()) + unR * n_if.y();
          uRg_for_L.z() = (uL.z() - unL * n_if.z()) + unR * n_if.z();

          const double rhoRg_L = rhoL;
          const double pRg_L = pR;
          double eRg_L{};
          double cRg_L{};
          vec<3> ruRg_L{};
          double reRg_L{};

          compute_state_from_rho_u_p<3>(eos_mm,
            mL,
            rhoRg_L,
            uRg_for_L,
            pRg_L,
            eRg_L,
            cRg_L,
            ruRg_L,
            reRg_L);

          const double f_r_T_L{ruL.y()};
          const double f_r_H_L{ruRg_L.y()};
          const vec<3> f_ru_T_L{
            ruL.y() * uL.x(), ruL.y() * uL.y() + pL, ruL.y() * uL.z()};
          const vec<3> f_ru_H_L{ruRg_L.y() * uRg_for_L.x(),
            ruRg_L.y() * uRg_for_L.y() + pRg_L,
            ruRg_L.y() * uRg_for_L.z()};
          const double f_rE_T_L{(reL + pL) * uL.y()};
          const double f_rE_H_L{(reRg_L + pRg_L) * uRg_for_L.y()};

          const double r_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_r_T_L,
              rhoRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_r_H_L,
              numerical_algorithms::field_tag::rho);

          const vec<3> ru_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_ru_T_L,
              ruRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_ru_H_L,
              numerical_algorithms::field_tag::ru);

          const double re_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_rE_T_L,
              reRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_rE_H_L,
              numerical_algorithms::field_tag::re);

          vec<3> n_rev{};
          n_rev.x() = -n_if.x();
          n_rev.y() = -n_if.y();
          n_rev.z() = -n_if.z();

          double unR_rev =
            uR.x() * n_rev.x() + uR.y() * n_rev.y() + uR.z() * n_rev.z();
          double unL_rev =
            uL.x() * n_rev.x() + uL.y() * n_rev.y() + uL.z() * n_rev.z();

          vec<3> uLg_for_R{};
          uLg_for_R.x() = (uR.x() - unR_rev * n_rev.x()) + unL_rev * n_rev.x();
          uLg_for_R.y() = (uR.y() - unR_rev * n_rev.y()) + unL_rev * n_rev.y();
          uLg_for_R.z() = (uR.z() - unR_rev * n_rev.z()) + unL_rev * n_rev.z();

          const double rhoLg_R = rhoR;
          const double pLg_R = pL;
          double eLg_R{};
          double cLg_R{};
          vec<3> ruLg_R{};
          double reLg_R{};

          compute_state_from_rho_u_p<3>(eos_mm,
            mR,
            rhoLg_R,
            uLg_for_R,
            pLg_R,
            eLg_R,
            cLg_R,
            ruLg_R,
            reLg_R);

          const double f_r_T_R{ruLg_R.y()};
          const double f_r_H_R{ruR.y()};
          const vec<3> f_ru_T_R{ruLg_R.y() * uLg_for_R.x(),
            ruLg_R.y() * uLg_for_R.y() + pLg_R,
            ruLg_R.y() * uLg_for_R.z()};
          const vec<3> f_ru_H_R{
            ruR.y() * uR.x(), ruR.y() * uR.y() + pR, ruR.y() * uR.z()};
          const double f_rE_T_R{(reLg_R + pLg_R) * uLg_for_R.y()};
          const double f_rE_H_R{(reR + pR) * uR.y()};

          const double r_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_r_T_R,
              rhoR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_r_H_R,
              numerical_algorithms::field_tag::rho);

          const vec<3> ru_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_ru_T_R,
              ruR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_ru_H_R,
              numerical_algorithms::field_tag::ru);

          const double re_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_rE_T_R,
              reR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_rE_H_R,
              numerical_algorithms::field_tag::re);

          r_f(i, j, k) = r_f_L;
          ru_f(i, j, k) = ru_f_L;
          re_f(i, j, k) = re_f_L;

          dt_mass_density(i, j - 1, k) -= inv_dx * r_f_L;
          dt_momentum_density(i, j - 1, k) -= inv_dx * ru_f_L;
          dt_total_energy_density(i, j - 1, k) -= inv_dx * re_f_L;

          dt_mass_density(i, j, k) += inv_dx * r_f_R;
          dt_momentum_density(i, j, k) += inv_dx * ru_f_R;
          dt_total_energy_density(i, j, k) += inv_dx * re_f_R;
        }
      };
    }
    else {
      const auto mdpolicy_cqq = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::corrector>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      s.executor().forall(kji, mdpolicy_cqq) {
        auto [k, j, i] = kji;
        const auto inv_dx = one_over_dx_i[2];

        const int mL = mat_id(i, j, k - 1);
        const int mR = mat_id(i, j, k);

        const double pLeft_wave =
          p_left(i, j, k - 1) +
          0.5 * r_left(i, j, k - 1) * g.z() / one_over_dx_i[2];
        const double pRight_wave =
          p_right(i, j, k) - 0.5 * r_right(i, j, k) * g.z() / one_over_dx_i[2];

        if(mL == mR) {
          const double f_r_T{ru_left(i, j, k - 1).z()};
          const double f_r_H{ru_right(i, j, k).z()};
          const vec<3> f_ru_T{
            ru_left(i, j, k - 1).z() * u_left(i, j, k - 1).x(),
            ru_left(i, j, k - 1).z() * u_left(i, j, k - 1).y(),
            ru_left(i, j, k - 1).z() * u_left(i, j, k - 1).z() + pLeft_wave};
          const vec<3> f_ru_H{ru_right(i, j, k).z() * u_right(i, j, k).x(),
            ru_right(i, j, k).z() * u_right(i, j, k).y(),
            ru_right(i, j, k).z() * u_right(i, j, k).z() + pRight_wave};
          const double f_rE_T{
            (re_left(i, j, k - 1) + pLeft_wave) * u_left(i, j, k - 1).z()};
          const double f_rE_H{
            (re_right(i, j, k) + pRight_wave) * u_right(i, j, k).z()};

          r_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              r_left(i, j, k - 1),
              r_left(i, j, k - 1),
              u_left(i, j, k - 1),
              re_left(i, j, k - 1),
              pLeft_wave,
              c_left(i, j, k - 1),
              f_r_T,
              r_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_r_H,
              numerical_algorithms::field_tag::rho);

          ru_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ru_left(i, j, k - 1),
              r_left(i, j, k - 1),
              u_left(i, j, k - 1),
              re_left(i, j, k - 1),
              pLeft_wave,
              c_left(i, j, k - 1),
              f_ru_T,
              ru_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_ru_H,
              numerical_algorithms::field_tag::ru);

          re_f(i, j, k) =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              re_left(i, j, k - 1),
              r_left(i, j, k - 1),
              u_left(i, j, k - 1),
              re_left(i, j, k - 1),
              pLeft_wave,
              c_left(i, j, k - 1),
              f_rE_T,
              re_right(i, j, k),
              r_right(i, j, k),
              u_right(i, j, k),
              re_right(i, j, k),
              pRight_wave,
              c_right(i, j, k),
              f_rE_H,
              numerical_algorithms::field_tag::re);

          dt_mass_density(i, j, k - 1) -= inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k - 1) -= inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k - 1) -= inv_dx * re_f(i, j, k);

          dt_mass_density(i, j, k) += inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k) += inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k) += inv_dx * re_f(i, j, k);
        }
        else {

          const double rhoL = std::max(r_left(i, j, k - 1), 1e-14);
          const vec<3> uL = u_left(i, j, k - 1);
          const double pL = std::max(pLeft_wave, 1e-14);
          const double cL = std::max(c_left(i, j, k - 1), 1e-14);
          const vec<3> ruL = ru_left(i, j, k - 1);
          const double reL = std::max(re_left(i, j, k - 1), 1e-14);

          const double rhoR = std::max(r_right(i, j, k), 1e-14);
          const vec<3> uR = u_right(i, j, k);
          const double pR = std::max(pRight_wave, 1e-14);
          const double cR = std::max(c_right(i, j, k), 1e-14);
          const vec<3> ruR = ru_right(i, j, k);
          const double reR = std::max(re_right(i, j, k), 1e-14);

          vec<3> n_if{};
          n_if.x() = 0.5 * (n_left(i, j, k - 1).x() + n_right(i, j, k).x());
          n_if.y() = 0.5 * (n_left(i, j, k - 1).y() + n_right(i, j, k).y());
          n_if.z() = 0.5 * (n_left(i, j, k - 1).z() + n_right(i, j, k).z());
          {
            const double nn = std::sqrt(
              n_if.x() * n_if.x() + n_if.y() * n_if.y() + n_if.z() * n_if.z());
            if(nn > 1e-14) {
              n_if.x() /= nn;
              n_if.y() /= nn;
              n_if.z() /= nn;
            }
            else {
              n_if.x() = 0.0;
              n_if.y() = 0.0;
              n_if.y() = 1.0;
            }
          }
          double unL =
            uL.x() * n_if.x() + uL.y() * n_if.y() + uL.z() * n_if.z();
          double unR =
            uR.x() * n_if.x() + uR.y() * n_if.y() + uR.z() * n_if.z();

          vec<3> uRg_for_L{};
          uRg_for_L.x() = (uL.x() - unL * n_if.x()) + unR * n_if.x();
          uRg_for_L.y() = (uL.y() - unL * n_if.y()) + unR * n_if.y();
          uRg_for_L.z() = (uL.z() - unL * n_if.z()) + unR * n_if.z();

          const double rhoRg_L = rhoL;
          const double pRg_L = pR;
          double eRg_L{};
          double cRg_L{};
          vec<3> ruRg_L{};
          double reRg_L{};

          compute_state_from_rho_u_p<3>(eos_mm,
            mL,
            rhoRg_L,
            uRg_for_L,
            pRg_L,
            eRg_L,
            cRg_L,
            ruRg_L,
            reRg_L);

          const double f_r_T_L{ruL.z()};
          const double f_r_H_L{ruRg_L.z()};
          const vec<3> f_ru_T_L{
            ruL.z() * uL.x(), ruL.z() * uL.y(), ruL.z() * uL.z() + pL};
          const vec<3> f_ru_H_L{ruRg_L.z() * uRg_for_L.x(),
            ruRg_L.z() * uRg_for_L.y(),
            ruRg_L.z() * uRg_for_L.z() + pRg_L};
          const double f_rE_T_L{(reL + pL) * uL.z()};
          const double f_rE_H_L{(reRg_L + pRg_L) * uRg_for_L.z()};

          const double r_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_r_T_L,
              rhoRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_r_H_L,
              numerical_algorithms::field_tag::rho);

          const vec<3> ru_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_ru_T_L,
              ruRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_ru_H_L,
              numerical_algorithms::field_tag::ru);

          const double re_f_L =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reL,
              rhoL,
              uL,
              reL,
              pL,
              cL,
              f_rE_T_L,
              reRg_L,
              rhoRg_L,
              uRg_for_L,
              reRg_L,
              pRg_L,
              cRg_L,
              f_rE_H_L,
              numerical_algorithms::field_tag::re);

          vec<3> n_rev{};
          n_rev.x() = -n_if.x();
          n_rev.y() = -n_if.y();
          n_rev.z() = -n_if.z();

          double unR_rev =
            uR.x() * n_rev.x() + uR.y() * n_rev.y() + uR.z() * n_rev.z();
          double unL_rev =
            uL.x() * n_rev.x() + uL.y() * n_rev.y() + uL.z() * n_rev.z();

          vec<3> uLg_for_R{};
          uLg_for_R.x() = (uR.x() - unR_rev * n_rev.x()) + unL_rev * n_rev.x();
          uLg_for_R.y() = (uR.y() - unR_rev * n_rev.y()) + unL_rev * n_rev.y();
          uLg_for_R.z() = (uR.z() - unR_rev * n_rev.z()) + unL_rev * n_rev.z();

          const double rhoLg_R = rhoR;
          const double pLg_R = pL;
          double eLg_R{};
          double cLg_R{};
          vec<3> ruLg_R{};
          double reLg_R{};

          compute_state_from_rho_u_p<3>(eos_mm,
            mR,
            rhoLg_R,
            uLg_for_R,
            pLg_R,
            eLg_R,
            cLg_R,
            ruLg_R,
            reLg_R);

          const double f_r_T_R{ruLg_R.z()};
          const double f_r_H_R{ruR.z()};
          const vec<3> f_ru_T_R{ruLg_R.z() * uLg_for_R.x(),
            ruLg_R.z() * uLg_for_R.y(),
            ruLg_R.z() * uLg_for_R.z() + pLg_R};
          const vec<3> f_ru_H_R{
            ruR.z() * uR.x(), ruR.z() * uR.y(), ruR.z() * uR.z() + pR};
          const double f_rE_T_R{(reLg_R + pLg_R) * uLg_for_R.z()};
          const double f_rE_H_R{(reR + pR) * uR.z()};

          const double r_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              rhoLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_r_T_R,
              rhoR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_r_H_R,
              numerical_algorithms::field_tag::rho);

          const vec<3> ru_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
              ruLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_ru_T_R,
              ruR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_ru_H_R,
              numerical_algorithms::field_tag::ru);

          const double re_f_R =
            numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
              reLg_R,
              rhoLg_R,
              uLg_for_R,
              reLg_R,
              pLg_R,
              cLg_R,
              f_rE_T_R,
              reR,
              rhoR,
              uR,
              reR,
              pR,
              cR,
              f_rE_H_R,
              numerical_algorithms::field_tag::re);

          r_f(i, j, k) = r_f_L;
          ru_f(i, j, k) = ru_f_L;
          re_f(i, j, k) = re_f_L;

          dt_mass_density(i, j, k - 1) -= inv_dx * r_f_L;
          dt_momentum_density(i, j, k - 1) -= inv_dx * ru_f_L;
          dt_total_energy_density(i, j, k - 1) -= inv_dx * re_f_L;

          dt_mass_density(i, j, k) += inv_dx * r_f_R;
          dt_momentum_density(i, j, k) += inv_dx * ru_f_R;
          dt_total_energy_density(i, j, k) += inv_dx * re_f_R;
        }
      };
    }
  }
}

} // namespace hard::tasks::gfm

#endif // HARD_MODULE_GFM_INTERFACE_FLUXES_HH
