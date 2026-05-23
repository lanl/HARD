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

// Helper structures for flux computation
template<std::size_t Dim>
struct flux_result {
  double r_flux;
  vec<Dim> ru_flux;
  double re_flux;
};

template<std::size_t Dim>
struct ghost_state {
  double rho;
  vec<Dim> u;
  double p;
  double e;
  double c;
  vec<Dim> ru;
  double re;
};

// Create artificial state for GFM cut-face.
// We need to provide complete and EOS-consistent state
// for the ghost side to HLLC (sound speed, total E den, mom).
template<std::size_t Dim>
void
compute_state_from_rho_u_p(const std::vector<eos::eos_wrapper> & eos_mm,
  int mat_id,
  ghost_state<Dim> & ghost) noexcept {

  constexpr double floor = 1e-16;

  flog_assert(std::isfinite(ghost.rho) && ghost.rho > floor,
    "compute_state_from_rho_u_p: non-finite/invalid ghost rho"
      << " rho=" << ghost.rho << " p=" << ghost.p << " mat_id=" << mat_id);

  flog_assert(std::isfinite(ghost.p) && ghost.p > floor,
    "compute_state_from_rho_u_p: non-finite/invalid ghost pressure"
      << " rho=" << ghost.rho << " p=" << ghost.p << " mat_id=" << mat_id);

  flog_assert(mat_id >= 0,
    "compute_state_from_rho_u_p: negative material id"
      << " mat_id=" << mat_id << " rho=" << ghost.rho << " p=" << ghost.p);

  const auto mat = static_cast<std::size_t>(mat_id);

  flog_assert(mat < eos_mm.size(),
    "compute_state_from_rho_u_p: material id out of EOS range"
      << " mat_id=" << mat_id << " eos_mm.size=" << eos_mm.size()
      << " rho=" << ghost.rho << " p=" << ghost.p);

  const auto & eos = eos_mm[mat];

  // Getting internal energy.
  ghost.e = hard::tasks::util::find_sie(eos, ghost.rho, ghost.p);

  flog_assert(std::isfinite(ghost.e) && ghost.e > 0.0,
    "compute_state_from_rho_u_p: invalid ghost sie"
      << " e=" << ghost.e << " rho=" << ghost.rho << " p=" << ghost.p
      << " mat_id=" << mat_id);

  // Compute sound speed.
  ghost.c = eos.cRhoSie(ghost.rho, ghost.e);

  flog_assert(std::isfinite(ghost.c) && ghost.c > 0.0,
    "compute_state_from_rho_u_p: invalid ghost sound speed"
      << " c=" << ghost.c << " e=" << ghost.e << " rho=" << ghost.rho
      << " p=" << ghost.p << " mat_id=" << mat_id);

  ghost.ru = ghost.rho * ghost.u;

  // Compute total energy density.
  ghost.re = ghost.rho * ghost.e + 0.5 * ghost.rho * ghost.u.norm_squared();

  flog_assert(std::isfinite(ghost.re) && ghost.re > 0.0,
    "compute_state_from_rho_u_p: invalid ghost total energy density"
      << " re=" << ghost.re << " rho=" << ghost.rho << " e=" << ghost.e
      << " p=" << ghost.p << " mat_id=" << mat_id);
}

// Helper function: compute all three HLLC fluxes at once
template<std::size_t Dim>
flux_result<Dim>
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
  double reR) noexcept {
  flux_result<Dim> result;

  // Compute flux terms
  const double f_r_T = ruL[face_axis];
  const double f_r_H = ruR[face_axis];

  vec<Dim> f_ru_T = ruL[face_axis] * uL;
  vec<Dim> f_ru_H = ruR[face_axis] * uR;
  f_ru_T[face_axis] += pL;
  f_ru_H[face_axis] += pR;

  const double f_rE_T = (reL + pL) * uL[face_axis];
  const double f_rE_H = (reR + pR) * uR[face_axis];

  // clang-format off
  result.r_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
      rhoL, rhoL, uL, reL, pL, cL, f_r_T, rhoR, rhoR, uR, reR,
      pR, cR, f_r_H, numerical_algorithms::field_tag::rho);
  // clang-format on

  // clang-format off
  result.ru_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
      ruL, rhoL, uL, reL, pL, cL, f_ru_T, ruR, rhoR, uR, reR, pR,
      cR, f_ru_H, numerical_algorithms::field_tag::ru);
  // clang-format on

  // clang-format off
  result.re_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
      reL, rhoL, uL, reL, pL, cL, f_rE_T, reR, rhoR, uR, reR,
      pR, cR, f_rE_H, numerical_algorithms::field_tag::re);
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

// Helper function: Validate whether extrapolated value usable
// as a ghost state for EOS/HLLC
template<std::size_t Dim>
bool
valid_extp_q_for_gfm(const hard::tasks::util::extp_q<Dim> & q) noexcept {
  return std::isfinite(q.rho()) && q.rho() > 0.0 && std::isfinite(q.p()) &&
         q.p() > 0.0 && std::isfinite(q.u().norm());
}

// Helper function: Check mat ids for valid EOS indices
inline bool
valid_material_id(int mat_id,
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {
  return mat_id >= 0 && static_cast<std::size_t>(mat_id) < eos_mm.size();
}

// Helper for raw ghost checks
template<std::size_t Dim>
void
assert_raw_ghost_state(const char * label,
  std::size_t face_axis,
  int mL,
  int mR,
  int mat_id,
  const ghost_state<Dim> & ghost) noexcept {

  flog_assert(std::isfinite(ghost.rho) && ghost.rho > 0.0,
    label << ": invalid raw ghost rho before EOS"
          << " face_axis=" << face_axis << " mL=" << mL << " mR=" << mR
          << " mat_id=" << mat_id << " rho=" << ghost.rho << " p=" << ghost.p
          << " e=" << ghost.e << " u_norm=" << ghost.u.norm());

  flog_assert(std::isfinite(ghost.p) && ghost.p > 0.0,
    label << ": invalid raw ghost pressure before EOS"
          << " face_axis=" << face_axis << " mL=" << mL << " mR=" << mR
          << " mat_id=" << mat_id << " rho=" << ghost.rho << " p=" << ghost.p
          << " e=" << ghost.e << " u_norm=" << ghost.u.norm());

  flog_assert(std::isfinite(ghost.u.norm()),
    label << ": invalid raw ghost velocity before EOS"
          << " face_axis=" << face_axis << " mL=" << mL << " mR=" << mR
          << " mat_id=" << mat_id << " rho=" << ghost.rho << " p=" << ghost.p
          << " e=" << ghost.e << " u_norm=" << ghost.u.norm());
}

// Helper function: compute ghost state for GFM cut face
template<std::size_t Dim>
ghost_state<Dim>
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

  ghost_state<Dim> ghost;
  ghost.rho = rho_owner;
  ghost.p = p_neighbor;
  ghost.u = u_ghost;

  compute_state_from_rho_u_p<Dim>(eos_mm, mat_id, ghost);

  return ghost;
}

// Helper function: Unpack ghost state and fill initial value from
// extrapolation.
template<std::size_t Dim>
ghost_state<Dim>
unpack_raw_ghost_state(const hard::tasks::util::extp_q<Dim> & q) noexcept {
  ghost_state<Dim> ghost;
  ghost.rho = q.rho();
  ghost.u = q.u();
  ghost.p = q.p();
  ghost.e = q.e();
  ghost.c = 0.0;
  ghost.ru = vec<Dim>{};
  ghost.re = 0.0;
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
  flux_result<Dim> & fluxes_L_out,
  flux_result<Dim> & fluxes_R_out,
  bool & is_cut_face) noexcept {

  is_cut_face = (mL != mR);

  if(!is_cut_face) {
    // Same material face
    auto fluxes = compute_all_HLLC_fluxes<Dim>(
      face_axis, rhoL, uL, pL, cL, ruL, reL, rhoR, uR, pR, cR, ruR, reR);

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
    auto ghostR_L = unpack_raw_ghost_state<Dim>(qGhostR_L);
    assert_raw_ghost_state<Dim>("ghostR_L", face_axis, mL, mR, mL, ghostR_L);
    compute_state_from_rho_u_p<Dim>(eos_mm, mL, ghostR_L);

    // Ghost state for the right owner problem:
    // right-material state extrapolated into the left side against right real
    // state.
    auto ghostL_R = unpack_raw_ghost_state<Dim>(qGhostL_R);
    assert_raw_ghost_state<Dim>("ghostL_R", face_axis, mL, mR, mR, ghostL_R);
    compute_state_from_rho_u_p<Dim>(eos_mm, mR, ghostL_R);

    // Left-owner flux
    // clang-format off
    fluxes_L_out = compute_all_HLLC_fluxes<Dim>(face_axis,
      rhoL, uL, pL, cL, ruL, reL, ghostR_L.rho, ghostR_L.u, ghostR_L.p,
      ghostR_L.c, ghostR_L.ru, ghostR_L.re);
    // Right-owner flux
    fluxes_R_out = compute_all_HLLC_fluxes<Dim>(face_axis,
      ghostL_R.rho, ghostL_R.u, ghostL_R.p, ghostL_R.c, ghostL_R.ru,
      ghostL_R.re, rhoR, uR, pR, cR, ruR, reR);
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

  // GFM switch
  constexpr bool enable_gfm = true;

  if constexpr(Dim == 1) {
    flog_assert(
      face_axis == 0, "1D executable can only compute fluxes along X axis");

    const auto inv_dx = one_over_dx_i[0];
    const auto dx = 1.0 / inv_dx;
    const double g_comp = g[0];

    // -------------------------------------------------------------------------
    // Compute face fluxes on dm::corrector
    //
    // For normal faces:
    //   r_f(i), ru_f(i), re_f(i) = regular HLLC flux.
    //
    // For GFM cut faces:
    //   r_f(i), ru_f(i), re_f(i) = left-owner GFM flux.
    //   Then add a correction to the right cell so it receives right-owner
    //   flux.
    // -------------------------------------------------------------------------
    for(auto i : m.template cells<ax::x, dm::corrector>()) {
      const int mL_raw = mat_id(i - 1);
      const int mR_raw = mat_id(i);

      const double pLeft_wave =
        p_left(i - 1) + 0.5 * r_left(i - 1) * g_comp * dx;
      const double pRight_wave = p_right(i) - 0.5 * r_right(i) * g_comp * dx;

      const double rhoL = r_left(i - 1);
      const double rhoR = r_right(i);
      const double pL = pLeft_wave;
      const double pR = pRight_wave;
      const double cL = c_left(i - 1);
      const double cR = c_right(i);
      const double reL = re_left(i - 1);
      const double reR = re_right(i);

      // check boundary face with valid materail
      const bool valid_material_pair =
        valid_material_id(mL_raw, eos_mm) && valid_material_id(mR_raw, eos_mm);

      // check boundary face with unfilled recon states
      const bool face_state_valid =
        std::isfinite(rhoL) && rhoL > 0.0 && std::isfinite(rhoR) &&
        rhoR > 0.0 && std::isfinite(pL) && pL > 0.0 && std::isfinite(pR) &&
        pR > 0.0 && std::isfinite(cL) && cL > 0.0 && std::isfinite(cR) &&
        cR > 0.0 && std::isfinite(reL) && reL > 0.0 && std::isfinite(reR) &&
        reR > 0.0;

      flog_assert(std::isfinite(rhoL) && rhoL >= 0.0,
        "rhoL(" << i << ") invalid: " << rhoL << " mL=" << mL_raw
                << " mR=" << mR_raw);
      flog_assert(std::isfinite(rhoR) && rhoR >= 0.0,
        "rhoR(" << i << ") invalid: " << rhoR << " mL=" << mL_raw
                << " mR=" << mR_raw);
      flog_assert(std::isfinite(pL) && pL >= 0.0,
        "pL(" << i << ") invalid: " << pL << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(pR) && pR >= 0.0,
        "pR(" << i << ") invalid: " << pR << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(cL) && cL >= 0.0,
        "cL(" << i << ") invalid: " << cL << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(cR) && cR >= 0.0,
        "cR(" << i << ") invalid: " << cR << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(reL) && reL >= 0.0,
        "reL(" << i << ") invalid: " << reL << " mL=" << mL_raw
               << " mR=" << mR_raw);
      flog_assert(std::isfinite(reR) && reR >= 0.0,
        "reR(" << i << ") invalid: " << reR << " mL=" << mL_raw
               << " mR=" << mR_raw);

      const bool qext_valid_for_cut =
        valid_extp_q_for_gfm<Dim>(q_ext_plus(i)) &&
        valid_extp_q_for_gfm<Dim>(q_ext_minus(i - 1));

      const bool use_gfm_cut = enable_gfm && valid_material_pair &&
                               face_state_valid && (mL_raw != mR_raw) &&
                               qext_valid_for_cut;

      // Force normal HLLC branch if not using GFM.
      // Preserve local left material ID for right material
      // if material 0 doesn't exist.
      const int mL_eff = mL_raw;
      const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      // const int mL_eff = mL_raw;
      // const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      flux_result<Dim> fluxes_L, fluxes_R;
      bool is_cut;

      // clang-format off
      process_single_face<Dim>(0,
        mL_eff, mR_eff, rhoL, u_left(i - 1), pL, cL, ru_left(i - 1), reL, rhoR,
        u_right(i), pR, cR, ru_right(i), reR, n_left(i - 1), n_right(i),
        q_ext_plus(i), q_ext_minus(i-1), eos_mm, r_f(i), ru_f(i),
        re_f(i), fluxes_L, fluxes_R, is_cut);
      // clang-format on

      // For a true GFM cut:
      //
      //   r_f/ru_f/re_f store the left-owner flux.
      //
      // The regular divergence update below gives:
      //   left cell  -= flux_L / dx
      //   right cell += flux_L / dx
      //
      // But the right cell should receive flux_R, so add:
      //   +(flux_R - flux_L) / dx
      //
      if(use_gfm_cut) {
        dt_mass_density(i) += inv_dx * (fluxes_R.r_flux - fluxes_L.r_flux);
        dt_momentum_density(i) +=
          inv_dx * (fluxes_R.ru_flux - fluxes_L.ru_flux);
        dt_total_energy_density(i) +=
          inv_dx * (fluxes_R.re_flux - fluxes_L.re_flux);
      }
    }

    // Apply flux divergence only on dm::quantities (like in regular hydro).
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      dt_mass_density(i) += inv_dx * (r_f(i) - r_f(i + 1));
      dt_momentum_density(i) += inv_dx * (ru_f(i) - ru_f(i + 1));
      dt_total_energy_density(i) += inv_dx * (re_f(i) - re_f(i + 1));
    }
  }
  // ---------------------------------------------------------------------------
  // 2D
  // ---------------------------------------------------------------------------
  else if constexpr(Dim == 2) {
    flog_assert(
      face_axis <= 1, "2D executable cannot compute fluxes along Z axis");

    // Index offsets for left cell: [face_axis][dimension]
    const std::array<std::array<int, 2>, 2> left_offset = {{
      {-1, 0}, // face_axis == 0: i-1, j
      {0, -1} // face_axis == 1: i, j-1
    }};

    const auto & offset = left_offset[face_axis];
    const auto inv_dx = one_over_dx_i[face_axis];
    const double g_comp = g[face_axis];

    auto mdpolicy = (face_axis == 0)
                      ? get_mdiota_policy(r_f,
                          m.template cells<ax::y, dm::quantities>(),
                          m.template cells<ax::x, dm::corrector>())
                      : get_mdiota_policy(r_f,
                          m.template cells<ax::y, dm::corrector>(),
                          m.template cells<ax::x, dm::quantities>());

    for(auto ji : mdpolicy) {
      auto [j, i] = ji;
      const int iL = i + offset[0];
      const int jL = j + offset[1];

      const int mL_raw = mat_id(iL, jL);
      const int mR_raw = mat_id(i, j);

      const double pLeft_wave =
        p_left(iL, jL) + 0.5 * r_left(iL, jL) * g_comp / inv_dx;
      const double pRight_wave =
        p_right(i, j) - 0.5 * r_right(i, j) * g_comp / inv_dx;

      const double rhoL = r_left(iL, jL);
      const double rhoR = r_right(i, j);
      const double pL = pLeft_wave;
      const double pR = pRight_wave;
      const double cL = c_left(iL, jL);
      const double cR = c_right(i, j);
      const double reL = re_left(iL, jL);
      const double reR = re_right(i, j);

      flog_assert(std::isfinite(rhoL) && rhoL >= 0.0,
        "rhoL(" << i << "," << j << ") invalid: " << rhoL << " mL=" << mL_raw
                << " mR=" << mR_raw);
      flog_assert(std::isfinite(rhoR) && rhoR >= 0.0,
        "rhoR(" << i << "," << j << ") invalid: " << rhoR << " mL=" << mL_raw
                << " mR=" << mR_raw);
      flog_assert(std::isfinite(pL) && pL >= 0.0,
        "pL(" << i << "," << j << ") invalid: " << pL << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(pR) && pR >= 0.0,
        "pR(" << i << "," << j << ") invalid: " << pR << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(cL) && cL >= 0.0,
        "cL(" << i << "," << j << ") invalid: " << cL << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(cR) && cR >= 0.0,
        "cR(" << i << "," << j << ") invalid: " << cR << " mL=" << mL_raw
              << " mR=" << mR_raw);
      flog_assert(std::isfinite(reL) && reL >= 0.0,
        "reL(" << i << "," << j << ") invalid: " << reL << " mL=" << mL_raw
               << " mR=" << mR_raw);
      flog_assert(std::isfinite(reR) && reR >= 0.0,
        "reR(" << i << "," << j << ") invalid: " << reR << " mL=" << mL_raw
               << " mR=" << mR_raw);

      const bool valid_material_pair =
        valid_material_id(mL_raw, eos_mm) && valid_material_id(mR_raw, eos_mm);

      // Strict check only for deciding whether GFM is allowed.
      const bool face_state_valid_for_gfm =
        std::isfinite(rhoL) && rhoL > 0.0 && std::isfinite(rhoR) &&
        rhoR > 0.0 && std::isfinite(pL) && pL > 0.0 && std::isfinite(pR) &&
        pR > 0.0 && std::isfinite(cL) && cL > 0.0 && std::isfinite(cR) &&
        cR > 0.0 && std::isfinite(reL) && reL > 0.0 && std::isfinite(reR) &&
        reR > 0.0;

      const bool qext_valid_for_cut =
        valid_extp_q_for_gfm<Dim>(q_ext_plus(i, j)) &&
        valid_extp_q_for_gfm<Dim>(q_ext_minus(iL, jL));

      const bool use_gfm_cut = enable_gfm && valid_material_pair &&
                               face_state_valid_for_gfm && (mL_raw != mR_raw) &&
                               qext_valid_for_cut;

      // Force normal HLLC branch if not using GFM.
      const int mL_eff = mL_raw;
      const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      flux_result<Dim> fluxes_L, fluxes_R;
      bool is_cut = false;
      // clang-format off
      process_single_face<Dim>(face_axis,
        mL_eff, mR_eff, rhoL, u_left(iL, jL), pL, cL, ru_left(iL, jL),
        reL, rhoR, u_right(i, j), pR, cR, ru_right(i, j), reR,
        n_left(iL, jL), n_right(i, j),
        q_ext_plus(i,j), q_ext_minus(iL,jL), eos_mm, r_f(i, j),
        ru_f(i, j), re_f(i, j), fluxes_L, fluxes_R, is_cut);
      // clang-format on

      if(use_gfm_cut) {
        dt_mass_density(i, j) += inv_dx * (fluxes_R.r_flux - fluxes_L.r_flux);
        dt_momentum_density(i, j) +=
          inv_dx * (fluxes_R.ru_flux - fluxes_L.ru_flux);
        dt_total_energy_density(i, j) +=
          inv_dx * (fluxes_R.re_flux - fluxes_L.re_flux);
      }
    }
    // Apply regular hydro-style divergence only over dm::quantities.
    auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto ji : mdpolicy_qq) {
      auto [j, i] = ji;

      if(face_axis == 0) {
        dt_mass_density(i, j) += inv_dx * (r_f(i, j) - r_f(i + 1, j));
        dt_momentum_density(i, j) += inv_dx * (ru_f(i, j) - ru_f(i + 1, j));
        dt_total_energy_density(i, j) += inv_dx * (re_f(i, j) - re_f(i + 1, j));
      }
      else {
        dt_mass_density(i, j) += inv_dx * (r_f(i, j) - r_f(i, j + 1));
        dt_momentum_density(i, j) += inv_dx * (ru_f(i, j) - ru_f(i, j + 1));
        dt_total_energy_density(i, j) += inv_dx * (re_f(i, j) - re_f(i, j + 1));
      }
    }
  }

  // ---------------------------------------------------------------------------
  // 3D
  // ---------------------------------------------------------------------------
  else {
    // Index offsets for left cell: [face_axis][dimension]
    const std::array<std::array<int, 3>, 3> left_offset = {{
      {-1, 0, 0}, // face_axis == 0: i-1, j, k
      {0, -1, 0}, // face_axis == 1: i, j-1, k
      {0, 0, -1} // face_axis == 2: i, j, k-1
    }};

    const auto & offset = left_offset[face_axis];
    const auto inv_dx = one_over_dx_i[face_axis];
    const double g_comp = g[face_axis];

    auto mdpolicy =
      (face_axis == 0)   ? get_mdiota_policy(r_f,
                           m.template cells<ax::z, dm::quantities>(),
                           m.template cells<ax::y, dm::quantities>(),
                           m.template cells<ax::x, dm::corrector>())
      : (face_axis == 1) ? get_mdiota_policy(r_f,
                             m.template cells<ax::z, dm::quantities>(),
                             m.template cells<ax::y, dm::corrector>(),
                             m.template cells<ax::x, dm::quantities>())
                         : get_mdiota_policy(r_f,
                             m.template cells<ax::z, dm::corrector>(),
                             m.template cells<ax::y, dm::quantities>(),
                             m.template cells<ax::x, dm::quantities>());

    for(auto kji : mdpolicy) {
      auto [k, j, i] = kji;
      const int iL = i + offset[0];
      const int jL = j + offset[1];
      const int kL = k + offset[2];

      const int mL_raw = mat_id(iL, jL, kL);
      const int mR_raw = mat_id(i, j, k);

      const double pLeft_wave =
        p_left(iL, jL, kL) + 0.5 * r_left(iL, jL, kL) * g_comp / inv_dx;
      const double pRight_wave =
        p_right(i, j, k) - 0.5 * r_right(i, j, k) * g_comp / inv_dx;

      const double rhoL = r_left(iL, jL, kL);
      const double rhoR = r_right(i, j, k);
      const double pL = pLeft_wave;
      const double pR = pRight_wave;
      const double cL = c_left(iL, jL, kL);
      const double cR = c_right(i, j, k);
      const double reL = re_left(iL, jL, kL);
      const double reR = re_right(i, j, k);

      flog_assert(std::isfinite(rhoL) && rhoL >= 0.0,
        "rhoL(" << i << "," << j << "," << k << ") invalid: " << rhoL
                << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(rhoR) && rhoR >= 0.0,
        "rhoR(" << i << "," << j << "," << k << ") invalid: " << rhoR
                << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(pL) && pL >= 0.0,
        "pL(" << i << "," << j << "," << k << ") invalid: " << pL
              << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(pR) && pR >= 0.0,
        "pR(" << i << "," << j << "," << k << ") invalid: " << pR
              << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(cL) && cL >= 0.0,
        "cL(" << i << "," << j << "," << k << ") invalid: " << cL
              << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(cR) && cR >= 0.0,
        "cR(" << i << "," << j << "," << k << ") invalid: " << cR
              << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(reL) && reL >= 0.0,
        "reL(" << i << "," << j << "," << k << ") invalid: " << reL
               << " mL=" << mL_raw << " mR=" << mR_raw);
      flog_assert(std::isfinite(reR) && reR >= 0.0,
        "reR(" << i << "," << j << "," << k << ") invalid: " << reR
               << " mL=" << mL_raw << " mR=" << mR_raw);
      const bool valid_material_pair =
        valid_material_id(mL_raw, eos_mm) && valid_material_id(mR_raw, eos_mm);

      // Strict check only for deciding whether GFM is allowed.
      const bool face_state_valid_for_gfm =
        std::isfinite(rhoL) && rhoL > 0.0 && std::isfinite(rhoR) &&
        rhoR > 0.0 && std::isfinite(pL) && pL > 0.0 && std::isfinite(pR) &&
        pR > 0.0 && std::isfinite(cL) && cL > 0.0 && std::isfinite(cR) &&
        cR > 0.0 && std::isfinite(reL) && reL > 0.0 && std::isfinite(reR) &&
        reR > 0.0;

      const bool qext_valid_for_cut =
        valid_extp_q_for_gfm<Dim>(q_ext_plus(i, j, k)) &&
        valid_extp_q_for_gfm<Dim>(q_ext_minus(iL, jL, kL));

      const bool use_gfm_cut = enable_gfm && valid_material_pair &&
                               face_state_valid_for_gfm && (mL_raw != mR_raw) &&
                               qext_valid_for_cut;

      // Force normal HLLC branch if not using GFM.
      const int mL_eff = mL_raw;
      const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      flux_result<Dim> fluxes_L, fluxes_R;
      bool is_cut = false;
      // clang-format off
      process_single_face<Dim>(face_axis,
        mL_eff, mR_eff, rhoL, u_left(iL, jL, kL), pL, cL, ru_left(iL, jL, kL),
        reL, rhoR, u_right(i, j, k), pR, cR, ru_right(i, j, k), reR,
        n_left(iL, jL, kL), n_right(i, j, k),
        q_ext_plus(i,j,k), q_ext_minus(iL,jL,kL), eos_mm, r_f(i, j, k),
        ru_f(i, j, k), re_f(i, j, k), fluxes_L, fluxes_R, is_cut);
      // clang-format on
      // GFM correction to the right/plus cell only.
      if(use_gfm_cut) {
        dt_mass_density(i, j, k) +=
          inv_dx * (fluxes_R.r_flux - fluxes_L.r_flux);
        dt_momentum_density(i, j, k) +=
          inv_dx * (fluxes_R.ru_flux - fluxes_L.ru_flux);
        dt_total_energy_density(i, j, k) +=
          inv_dx * (fluxes_R.re_flux - fluxes_L.re_flux);
      }
    }

    // Apply regular hydro-style divergence only over dm::quantities.
    auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto kji : mdpolicy_qqq) {
      auto [k, j, i] = kji;

      if(face_axis == 0) {
        dt_mass_density(i, j, k) += inv_dx * (r_f(i, j, k) - r_f(i + 1, j, k));
        dt_momentum_density(i, j, k) +=
          inv_dx * (ru_f(i, j, k) - ru_f(i + 1, j, k));
        dt_total_energy_density(i, j, k) +=
          inv_dx * (re_f(i, j, k) - re_f(i + 1, j, k));
      }
      else if(face_axis == 1) {
        dt_mass_density(i, j, k) += inv_dx * (r_f(i, j, k) - r_f(i, j + 1, k));
        dt_momentum_density(i, j, k) +=
          inv_dx * (ru_f(i, j, k) - ru_f(i, j + 1, k));
        dt_total_energy_density(i, j, k) +=
          inv_dx * (re_f(i, j, k) - re_f(i, j + 1, k));
      }
      else {
        dt_mass_density(i, j, k) += inv_dx * (r_f(i, j, k) - r_f(i, j, k + 1));
        dt_momentum_density(i, j, k) +=
          inv_dx * (ru_f(i, j, k) - ru_f(i, j, k + 1));
        dt_total_energy_density(i, j, k) +=
          inv_dx * (re_f(i, j, k) - re_f(i, j, k + 1));
      }
    }
  }
}

} // namespace hard::tasks::gfm

#endif // HARD_MODULE_GFM_INTERFACE_FLUXES_HH
