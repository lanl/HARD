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

  ru = rho_f * u;
  // Compute total energy density
  rE = rho_f * e + 0.5 * rho_f * u.norm_squared();
}

// Helper structures for flux computation
template<std::size_t Dim>
struct FluxResult {
  double r_flux;
  vec<Dim> ru_flux;
  double re_flux;
};

template<std::size_t Dim>
struct GhostState {
  double rho;
  vec<Dim> u;
  double p;
  double e;
  double c;
  vec<Dim> ru;
  double re;
};

template<std::size_t Dim>
struct FluxTerms {
  double f_r_T, f_r_H;
  vec<Dim> f_ru_T, f_ru_H;
  double f_rE_T, f_rE_H;
};

// Helper function: compute flux terms for HLLC solver
template<std::size_t Dim>
FluxTerms<Dim>
compute_flux_terms(std::size_t face_axis,
  const vec<Dim> & ruL,
  const vec<Dim> & uL,
  double pL,
  double reL,
  const vec<Dim> & ruR,
  const vec<Dim> & uR,
  double pR,
  double reR) noexcept {
  FluxTerms<Dim> terms;

  terms.f_r_T = ruL[face_axis];
  terms.f_r_H = ruR[face_axis];

  terms.f_ru_T = ruL[face_axis] * uL;
  terms.f_ru_H = ruR[face_axis] * uR;
  terms.f_ru_T[face_axis] += pL;
  terms.f_ru_H[face_axis] += pR;

  terms.f_rE_T = (reL + pL) * uL[face_axis];
  terms.f_rE_H = (reR + pR) * uR[face_axis];

  return terms;
}

// Helper function: compute all three HLLC fluxes at once
template<std::size_t Dim>
FluxResult<Dim>
compute_all_HLLC_fluxes(std::size_t face_axis,
  double rhoL,
  const vec<Dim> & uL,
  double pL,
  double cL,
  const vec<Dim> & ruL,
  double reL,
  double rhoR,
  const vec<Dim> & uR,
  double pR,
  double cR,
  const vec<Dim> & ruR,
  double reR,
  const FluxTerms<Dim> & terms) noexcept {
  FluxResult<Dim> result;

  // clang-format off
  result.r_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
      rhoL, rhoL, uL, reL, pL, cL, terms.f_r_T, rhoR, rhoR, uR, reR,
      pR, cR, terms.f_r_H, numerical_algorithms::field_tag::rho);
  // clang-format on

  // clang-format off
  result.ru_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
      ruL, rhoL, uL, reL, pL, cL, terms.f_ru_T, ruR, rhoR, uR, reR, pR,
      cR, terms.f_ru_H, numerical_algorithms::field_tag::ru);
  // clang-format on

  // clang-format off
  result.re_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
      reL, rhoL, uL, reL, pL, cL, terms.f_rE_T, reR, rhoR, uR, reR,
      pR, cR, terms.f_rE_H, numerical_algorithms::field_tag::re);
  // clang-format on

  return result;
}

// Helper function: normalize interface normal from left/right views
template<std::size_t Dim>
vec<Dim>
normalize_interface_normal(const vec<Dim> & n_left,
  const vec<Dim> & n_right,
  std::size_t face_axis) noexcept {
  vec<Dim> n_if = 0.5 * (n_left + n_right);
  const double nn = std::sqrt(n_if.norm_squared());

  if(nn > 1e-14) {
    n_if = n_if * (1.0 / nn);
  }
  else {
    n_if = vec<Dim>{};
    n_if[face_axis] = 1.0;
  }

  return n_if;
}

// Helper function: compute ghost state for GFM cut face
template<std::size_t Dim>
GhostState<Dim>
compute_ghost_state(const std::vector<eos::eos_wrapper> & eos_mm,
  int mat_id,
  const vec<Dim> & n_if,
  double rho_owner,
  const vec<Dim> & u_owner,
  const vec<Dim> & u_neighbor,
  double p_neighbor) noexcept {

  double un_owner = 0.0;
  double un_neighbor = 0.0;
  for(std::size_t d = 0; d < Dim; ++d) {
    un_owner += u_owner[d] * n_if[d];
    un_neighbor += u_neighbor[d] * n_if[d];
  }

  const vec<Dim> u_ghost = u_owner - un_owner * n_if + un_neighbor * n_if;

  GhostState<Dim> ghost;
  ghost.rho = rho_owner;
  ghost.p = p_neighbor;
  ghost.u = u_ghost;

  // clang-format off
  compute_state_from_rho_u_p<Dim>(eos_mm,
    mat_id, ghost.rho, ghost.u, ghost.p, ghost.e, ghost.c, ghost.ru,
    ghost.re);
  // clang-format on

  return ghost;
}

// Helper function: build a complete EOS-consistent ghost state from
// an extrapolated packed primitive state.
template<std::size_t Dim>
GhostState<Dim>
compute_ghost_state_from_extp_q(const std::vector<eos::eos_wrapper> & eos_mm,
  int mat_id,
  const hard::tasks::util::extp_q<Dim> & q) noexcept {

  using ids = hard::tasks::util::extp_q_ids<Dim>;

  GhostState<Dim> ghost;

  ghost.rho = std::max(q[ids::rho], 1e-14);

  ghost.u = vec<Dim>(0.0);
  ghost.u.x() = q[ids::vel + 0];

  if constexpr(Dim >= 2) {
    ghost.u.y() = q[ids::vel + 1];
  }

  if constexpr(Dim >= 3) {
    ghost.u.z() = q[ids::vel + 2];
  }

  ghost.p = std::max(q[ids::pressure], 1e-14);

  // Recompute EOS-consistent state.
  compute_state_from_rho_u_p<Dim>(eos_mm,
    mat_id,
    ghost.rho,
    ghost.u,
    ghost.p,
    ghost.e,
    ghost.c,
    ghost.ru,
    ghost.re);

  return ghost;
}

// Helper function: process a single face (same material or GFM cut face)
template<std::size_t Dim>
void
process_single_face(std::size_t face_axis,
  int mL,
  int mR,
  double rhoL,
  const vec<Dim> & uL,
  double pL,
  double cL,
  const vec<Dim> & ruL,
  double reL,
  double rhoR,
  const vec<Dim> & uR,
  double pR,
  double cR,
  const vec<Dim> & ruR,
  double reR,
  const vec<Dim> & n_left,
  const vec<Dim> & n_right,
  const hard::tasks::util::extp_q<Dim> & qGhostR_L,
  const hard::tasks::util::extp_q<Dim> & qGhostL_R,
  const std::vector<eos::eos_wrapper> & eos_mm,
  double & r_flux_out,
  vec<Dim> & ru_flux_out,
  double & re_flux_out,
  FluxResult<Dim> & fluxes_L_out,
  FluxResult<Dim> & fluxes_R_out,
  bool & is_cut_face) noexcept {

  is_cut_face = (mL != mR);

  if(mL == mR) {
    // Same material face
    auto terms =
      compute_flux_terms<Dim>(face_axis, ruL, uL, pL, reL, ruR, uR, pR, reR);
    auto fluxes = compute_all_HLLC_fluxes<Dim>(
      face_axis, rhoL, uL, pL, cL, ruL, reL, rhoR, uR, pR, cR, ruR, reR, terms);

    r_flux_out = fluxes.r_flux;
    ru_flux_out = fluxes.ru_flux;
    re_flux_out = fluxes.re_flux;
  }
  else {
    // GFM cut face
    const vec<Dim> n_if =
      normalize_interface_normal(n_left, n_right, face_axis);

    // Ghost state for the left owner problem:
    // left real state against left-material state extrapolated into the right
    // side.
    auto ghostR_L = compute_ghost_state_from_extp_q<Dim>(eos_mm, mL, qGhostR_L);

    // Ghost state for the right owner problem:
    // right-material state extrapolated into the left side against right real
    // state.
    auto ghostL_R = compute_ghost_state_from_extp_q<Dim>(eos_mm, mR, qGhostL_R);

    // Left-owner flux
    // clang-format off
    auto terms_L = compute_flux_terms<Dim>(face_axis,
      ruL, uL, pL, reL, ghostR_L.ru, ghostR_L.u, ghostR_L.p, ghostR_L.re);
    fluxes_L_out = compute_all_HLLC_fluxes<Dim>(face_axis,
      rhoL, uL, pL, cL, ruL, reL, ghostR_L.rho, ghostR_L.u, ghostR_L.p,
      ghostR_L.c, ghostR_L.ru, ghostR_L.re, terms_L);
    // Right-owner flux
    auto terms_R = compute_flux_terms<Dim>(face_axis,
      ghostL_R.ru, ghostL_R.u, ghostL_R.p, ghostL_R.re, ruR, uR, pR, reR);
    fluxes_R_out = compute_all_HLLC_fluxes<Dim>(face_axis,
      ghostL_R.rho, ghostL_R.u, ghostL_R.p, ghostL_R.c, ghostL_R.ru,
      ghostL_R.re, rhoR, uR, pR, cR, ruR, reR, terms_R);
    // clang-format on

    // Store left-view flux
    r_flux_out = fluxes_L_out.r_flux;
    ru_flux_out = fluxes_L_out.ru_flux;
    re_flux_out = fluxes_L_out.re_flux;
  }
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
  typename field<hard::tasks::util::extp_q<Dim>>::template accessor<ro, ro>
    q_ext_plus_a,
  typename field<hard::tasks::util::extp_q<Dim>>::template accessor<ro, ro>
    q_ext_minus_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> r_f_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ru_f_a,
  field<double>::accessor<wo, ro> re_f_a,
  // time derivative
  typename RK<Dim>::template accessor<rw, na> rk_dt_a,
  typename single<vec<Dim>>::template accessor<ro> g_acc,
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {

  auto g = *g_acc;

  auto [r_right, r_left] = faces<Dim>::mdcolex(m, r_face_a);
  auto [u_right, u_left] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [p_right, p_left] = faces<Dim>::mdcolex(m, pFace_a);
  auto [c_right, c_left] = faces<Dim>::mdcolex(m, cFace_a);
  auto [ru_right, ru_left] = faces_vec<Dim>::mdcolex(m, ru_face_a);
  auto [re_right, re_left] = faces<Dim>::mdcolex(m, re_face_a);
  auto [n_right, n_left] = faces_vec<Dim>::mdcolex(m, normal_f_a);

  auto mat_id = m.template mdcolex<is::cells>(material_a);
  auto q_ext_plus = m.template mdcolex<is::cells>(q_ext_plus_a);
  auto q_ext_minus = m.template mdcolex<is::cells>(q_ext_minus_a);

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

  // ---------------------------------------------------------------------------
  // 1D
  // ---------------------------------------------------------------------------
  if constexpr(Dim == 1) {
    flog_assert(
      face_axis == 0, "1D executable can only compute fluxes along X axis");

    for(auto i : m.template cells<ax::x, dm::corrector>()) {
      const auto dx = m.template delta<ax::x>();
      const auto inv_dx = one_over_dx_i[0];

      const int mL = mat_id(i - 1);
      const int mR = mat_id(i);

      const double pLeft_wave =
        p_left(i - 1) + 0.5 * r_left(i - 1) * g.x() * dx;
      const double pRight_wave = p_right(i) - 0.5 * r_right(i) * g.x() * dx;

      const double rhoL = std::max(r_left(i - 1), 1e-14);
      const double rhoR = std::max(r_right(i), 1e-14);
      const double pL = std::max(pLeft_wave, 1e-14);
      const double pR = std::max(pRight_wave, 1e-14);
      const double cL = std::max(c_left(i - 1), 1e-14);
      const double cR = std::max(c_right(i), 1e-14);
      const double reL = std::max(re_left(i - 1), 1e-14);
      const double reR = std::max(re_right(i), 1e-14);

      FluxResult<Dim> fluxes_L, fluxes_R;
      bool is_cut;
      // clang-format off
      process_single_face<Dim>(face_axis,
        mL, mR, rhoL, u_left(i - 1), pL, cL, ru_left(i - 1), reL, rhoR,
        u_right(i), pR, cR, ru_right(i), reR, n_left(i - 1), n_right(i),
	q_ext_plus(i), q_ext_minus(i-1),
        eos_mm, r_f(i), ru_f(i), re_f(i), fluxes_L, fluxes_R, is_cut);
      // clang-format on

      if(!is_cut) {
        dt_mass_density(i - 1) -= inv_dx * r_f(i);
        dt_momentum_density(i - 1) -= inv_dx * ru_f(i);
        dt_total_energy_density(i - 1) -= inv_dx * re_f(i);

        dt_mass_density(i) += inv_dx * r_f(i);
        dt_momentum_density(i) += inv_dx * ru_f(i);
        dt_total_energy_density(i) += inv_dx * re_f(i);
      }
      else {
        dt_mass_density(i - 1) -= inv_dx * fluxes_L.r_flux;
        dt_momentum_density(i - 1) -= inv_dx * fluxes_L.ru_flux;
        dt_total_energy_density(i - 1) -= inv_dx * fluxes_L.re_flux;

        dt_mass_density(i) += inv_dx * fluxes_R.r_flux;
        dt_momentum_density(i) += inv_dx * fluxes_R.ru_flux;
        dt_total_energy_density(i) += inv_dx * fluxes_R.re_flux;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // 2D
  // ---------------------------------------------------------------------------
  else if constexpr(Dim == 2) {
    flog_assert(
      face_axis <= 1, "2D executable cannot compute fluxes along Z axis");

    if(face_axis == 0) {
      auto mdpolicy_qc = get_mdiota_policy(r_f,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::corrector>());

      for(auto ji : mdpolicy_qc) {
        auto [j, i] = ji;
        const auto inv_dx = one_over_dx_i[0];

        const int mL = mat_id(i - 1, j);
        const int mR = mat_id(i, j);

        const double pLeft_wave =
          p_left(i - 1, j) + 0.5 * r_left(i - 1, j) * g.x() / one_over_dx_i[0];
        const double pRight_wave =
          p_right(i, j) - 0.5 * r_right(i, j) * g.x() / one_over_dx_i[0];

        const double rhoL = std::max(r_left(i - 1, j), 1e-14);
        const double rhoR = std::max(r_right(i, j), 1e-14);
        const double pL = std::max(pLeft_wave, 1e-14);
        const double pR = std::max(pRight_wave, 1e-14);
        const double cL = std::max(c_left(i - 1, j), 1e-14);
        const double cR = std::max(c_right(i, j), 1e-14);
        const double reL = std::max(re_left(i - 1, j), 1e-14);
        const double reR = std::max(re_right(i, j), 1e-14);

        FluxResult<Dim> fluxes_L, fluxes_R;
        bool is_cut;
        // clang-format off
        process_single_face<Dim>(face_axis,
          mL, mR, rhoL, u_left(i - 1, j), pL, cL, ru_left(i - 1, j),
          reL, rhoR, u_right(i, j), pR, cR, ru_right(i, j), reR,
          n_left(i - 1, j), n_right(i, j),  
	  q_ext_plus(i,j), q_ext_minus(i-1,j),eos_mm, r_f(i, j),
          ru_f(i, j), re_f(i, j), fluxes_L, fluxes_R, is_cut);
        // clang-format on

        if(!is_cut) {
          dt_mass_density(i - 1, j) -= inv_dx * r_f(i, j);
          dt_momentum_density(i - 1, j) -= inv_dx * ru_f(i, j);
          dt_total_energy_density(i - 1, j) -= inv_dx * re_f(i, j);

          dt_mass_density(i, j) += inv_dx * r_f(i, j);
          dt_momentum_density(i, j) += inv_dx * ru_f(i, j);
          dt_total_energy_density(i, j) += inv_dx * re_f(i, j);
        }
        else {
          dt_mass_density(i - 1, j) -= inv_dx * fluxes_L.r_flux;
          dt_momentum_density(i - 1, j) -= inv_dx * fluxes_L.ru_flux;
          dt_total_energy_density(i - 1, j) -= inv_dx * fluxes_L.re_flux;

          dt_mass_density(i, j) += inv_dx * fluxes_R.r_flux;
          dt_momentum_density(i, j) += inv_dx * fluxes_R.ru_flux;
          dt_total_energy_density(i, j) += inv_dx * fluxes_R.re_flux;
        }
      }
    }
    else if(face_axis == 1) {
      auto mdpolicy_cq = get_mdiota_policy(r_f,
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      for(auto ji : mdpolicy_cq) {
        auto [j, i] = ji;
        const auto inv_dx = one_over_dx_i[1];

        const int mL = mat_id(i, j - 1);
        const int mR = mat_id(i, j);

        const double pLeft_wave =
          p_left(i, j - 1) + 0.5 * r_left(i, j - 1) * g.y() / one_over_dx_i[1];
        const double pRight_wave =
          p_right(i, j) - 0.5 * r_right(i, j) * g.y() / one_over_dx_i[1];

        const double rhoL = std::max(r_left(i, j - 1), 1e-14);
        const double rhoR = std::max(r_right(i, j), 1e-14);
        const double pL = std::max(pLeft_wave, 1e-14);
        const double pR = std::max(pRight_wave, 1e-14);
        const double cL = std::max(c_left(i, j - 1), 1e-14);
        const double cR = std::max(c_right(i, j), 1e-14);
        const double reL = std::max(re_left(i, j - 1), 1e-14);
        const double reR = std::max(re_right(i, j), 1e-14);

        FluxResult<Dim> fluxes_L, fluxes_R;
        bool is_cut;
        // clang-format off
        process_single_face<Dim>(face_axis,
          mL, mR, rhoL, u_left(i, j - 1), pL, cL, ru_left(i, j - 1),
          reL, rhoR, u_right(i, j), pR, cR, ru_right(i, j), reR,
          n_left(i, j - 1), n_right(i, j),
	  q_ext_plus(i,j), q_ext_minus(i,j-1),eos_mm, r_f(i, j),
          ru_f(i, j), re_f(i, j), fluxes_L, fluxes_R, is_cut);
        // clang-format on

        if(!is_cut) {
          dt_mass_density(i, j - 1) -= inv_dx * r_f(i, j);
          dt_momentum_density(i, j - 1) -= inv_dx * ru_f(i, j);
          dt_total_energy_density(i, j - 1) -= inv_dx * re_f(i, j);

          dt_mass_density(i, j) += inv_dx * r_f(i, j);
          dt_momentum_density(i, j) += inv_dx * ru_f(i, j);
          dt_total_energy_density(i, j) += inv_dx * re_f(i, j);
        }
        else {
          dt_mass_density(i, j - 1) -= inv_dx * fluxes_L.r_flux;
          dt_momentum_density(i, j - 1) -= inv_dx * fluxes_L.ru_flux;
          dt_total_energy_density(i, j - 1) -= inv_dx * fluxes_L.re_flux;

          dt_mass_density(i, j) += inv_dx * fluxes_R.r_flux;
          dt_momentum_density(i, j) += inv_dx * fluxes_R.ru_flux;
          dt_total_energy_density(i, j) += inv_dx * fluxes_R.re_flux;
        }
      }
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

      for(auto kji : mdpolicy_qqc) {
        auto [k, j, i] = kji;
        const auto inv_dx = one_over_dx_i[0];

        const int mL = mat_id(i - 1, j, k);
        const int mR = mat_id(i, j, k);

        const double pLeft_wave =
          p_left(i - 1, j, k) +
          0.5 * r_left(i - 1, j, k) * g.x() / one_over_dx_i[0];
        const double pRight_wave =
          p_right(i, j, k) - 0.5 * r_right(i, j, k) * g.x() / one_over_dx_i[0];

        const double rhoL = std::max(r_left(i - 1, j, k), 1e-14);
        const double rhoR = std::max(r_right(i, j, k), 1e-14);
        const double pL = std::max(pLeft_wave, 1e-14);
        const double pR = std::max(pRight_wave, 1e-14);
        const double cL = std::max(c_left(i - 1, j, k), 1e-14);
        const double cR = std::max(c_right(i, j, k), 1e-14);
        const double reL = std::max(re_left(i - 1, j, k), 1e-14);
        const double reR = std::max(re_right(i, j, k), 1e-14);

        FluxResult<Dim> fluxes_L, fluxes_R;
        bool is_cut;
        // clang-format off
        process_single_face<Dim>(face_axis,
          mL, mR, rhoL, u_left(i - 1, j, k), pL, cL, ru_left(i - 1, j, k),
          reL, rhoR, u_right(i, j, k), pR, cR, ru_right(i, j, k), reR,
          n_left(i - 1, j, k), n_right(i, j, k), 
	  q_ext_plus(i,j,k), q_ext_minus(i-1,j,k),eos_mm, r_f(i, j, k),
          ru_f(i, j, k), re_f(i, j, k), fluxes_L, fluxes_R, is_cut);
        // clang-format on

        if(!is_cut) {
          dt_mass_density(i - 1, j, k) -= inv_dx * r_f(i, j, k);
          dt_momentum_density(i - 1, j, k) -= inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i - 1, j, k) -= inv_dx * re_f(i, j, k);

          dt_mass_density(i, j, k) += inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k) += inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k) += inv_dx * re_f(i, j, k);
        }
        else {
          dt_mass_density(i - 1, j, k) -= inv_dx * fluxes_L.r_flux;
          dt_momentum_density(i - 1, j, k) -= inv_dx * fluxes_L.ru_flux;
          dt_total_energy_density(i - 1, j, k) -= inv_dx * fluxes_L.re_flux;

          dt_mass_density(i, j, k) += inv_dx * fluxes_R.r_flux;
          dt_momentum_density(i, j, k) += inv_dx * fluxes_R.ru_flux;
          dt_total_energy_density(i, j, k) += inv_dx * fluxes_R.re_flux;
        }
      }
    }
    else if(face_axis == 1) {
      const auto mdpolicy_qcq = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::corrector>(),
        m.template cells<ax::x, dm::quantities>());

      for(auto kji : mdpolicy_qcq) {
        auto [k, j, i] = kji;
        const auto inv_dx = one_over_dx_i[1];

        const int mL = mat_id(i, j - 1, k);
        const int mR = mat_id(i, j, k);

        const double pLeft_wave =
          p_left(i, j - 1, k) +
          0.5 * r_left(i, j - 1, k) * g.y() / one_over_dx_i[1];
        const double pRight_wave =
          p_right(i, j, k) - 0.5 * r_right(i, j, k) * g.y() / one_over_dx_i[1];

        const double rhoL = std::max(r_left(i, j - 1, k), 1e-14);
        const double rhoR = std::max(r_right(i, j, k), 1e-14);
        const double pL = std::max(pLeft_wave, 1e-14);
        const double pR = std::max(pRight_wave, 1e-14);
        const double cL = std::max(c_left(i, j - 1, k), 1e-14);
        const double cR = std::max(c_right(i, j, k), 1e-14);
        const double reL = std::max(re_left(i, j - 1, k), 1e-14);
        const double reR = std::max(re_right(i, j, k), 1e-14);

        FluxResult<Dim> fluxes_L, fluxes_R;
        bool is_cut;
        // clang-format off
        process_single_face<Dim>(face_axis,
          mL, mR, rhoL, u_left(i, j - 1, k), pL, cL, ru_left(i, j - 1, k),
          reL, rhoR, u_right(i, j, k), pR, cR, ru_right(i, j, k), reR,
          n_left(i, j - 1, k), n_right(i, j, k), 
	  q_ext_plus(i,j,k), q_ext_minus(i,j-1,k),eos_mm, r_f(i, j, k),
          ru_f(i, j, k), re_f(i, j, k), fluxes_L, fluxes_R, is_cut);
        // clang-format on

        if(!is_cut) {
          dt_mass_density(i, j - 1, k) -= inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j - 1, k) -= inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j - 1, k) -= inv_dx * re_f(i, j, k);

          dt_mass_density(i, j, k) += inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k) += inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k) += inv_dx * re_f(i, j, k);
        }
        else {
          dt_mass_density(i, j - 1, k) -= inv_dx * fluxes_L.r_flux;
          dt_momentum_density(i, j - 1, k) -= inv_dx * fluxes_L.ru_flux;
          dt_total_energy_density(i, j - 1, k) -= inv_dx * fluxes_L.re_flux;

          dt_mass_density(i, j, k) += inv_dx * fluxes_R.r_flux;
          dt_momentum_density(i, j, k) += inv_dx * fluxes_R.ru_flux;
          dt_total_energy_density(i, j, k) += inv_dx * fluxes_R.re_flux;
        }
      }
    }
    else {
      const auto mdpolicy_cqq = get_mdiota_policy(r_f,
        m.template cells<ax::z, dm::corrector>(),
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());

      for(auto kji : mdpolicy_cqq) {
        auto [k, j, i] = kji;
        const auto inv_dx = one_over_dx_i[2];

        const int mL = mat_id(i, j, k - 1);
        const int mR = mat_id(i, j, k);

        const double pLeft_wave =
          p_left(i, j, k - 1) +
          0.5 * r_left(i, j, k - 1) * g.z() / one_over_dx_i[2];
        const double pRight_wave =
          p_right(i, j, k) - 0.5 * r_right(i, j, k) * g.z() / one_over_dx_i[2];

        const double rhoL = std::max(r_left(i, j, k - 1), 1e-14);
        const double rhoR = std::max(r_right(i, j, k), 1e-14);
        const double pL = std::max(pLeft_wave, 1e-14);
        const double pR = std::max(pRight_wave, 1e-14);
        const double cL = std::max(c_left(i, j, k - 1), 1e-14);
        const double cR = std::max(c_right(i, j, k), 1e-14);
        const double reL = std::max(re_left(i, j, k - 1), 1e-14);
        const double reR = std::max(re_right(i, j, k), 1e-14);

        FluxResult<Dim> fluxes_L, fluxes_R;
        bool is_cut;
        // clang-format off
        process_single_face<Dim>(face_axis,
          mL, mR, rhoL, u_left(i, j, k - 1), pL, cL, ru_left(i, j, k - 1),
          reL, rhoR, u_right(i, j, k), pR, cR, ru_right(i, j, k), reR,
          n_left(i, j, k - 1), n_right(i, j, k), 
	  q_ext_plus(i,j,k), q_ext_minus(i,j,k-1),eos_mm, r_f(i, j, k),
          ru_f(i, j, k), re_f(i, j, k), fluxes_L, fluxes_R, is_cut);
        // clang-format on

        if(!is_cut) {
          dt_mass_density(i, j, k - 1) -= inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k - 1) -= inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k - 1) -= inv_dx * re_f(i, j, k);

          dt_mass_density(i, j, k) += inv_dx * r_f(i, j, k);
          dt_momentum_density(i, j, k) += inv_dx * ru_f(i, j, k);
          dt_total_energy_density(i, j, k) += inv_dx * re_f(i, j, k);
        }
        else {
          dt_mass_density(i, j, k - 1) -= inv_dx * fluxes_L.r_flux;
          dt_momentum_density(i, j, k - 1) -= inv_dx * fluxes_L.ru_flux;
          dt_total_energy_density(i, j, k - 1) -= inv_dx * fluxes_L.re_flux;

          dt_mass_density(i, j, k) += inv_dx * fluxes_R.r_flux;
          dt_momentum_density(i, j, k) += inv_dx * fluxes_R.ru_flux;
          dt_total_energy_density(i, j, k) += inv_dx * fluxes_R.re_flux;
        }
      }
    }
  }
}

} // namespace hard::tasks::gfm

#endif // HARD_MODULE_GFM_INTERFACE_FLUXES_HH
