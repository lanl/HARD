#ifndef HARD_MODULE_GFM_INTERFACE_FLUXES_HH
#define HARD_MODULE_GFM_INTERFACE_FLUXES_HH

#include "../../hydro/numerical_algorithms/riemann_solvers.hh"
#include "../../rad/tasks/utils.hh"
#include "../modules/common/tasks/utils.hh"
#include "../modules/spec/eos.hh"
#include "../state.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>

namespace hard::tasks::gfm {
// GFM method toggle: Set to true for original gfm. Set to false for Riemann gfm
// -----------------------------------------------------------------------------

inline constexpr bool use_original_gfm = true;

// -----------------------------------------------------------------------------
// Helper structures for flux computation
// -----------------------------------------------------------------------------

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

template<std::size_t D>
void
diagnose_rk_rhs_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename RK<D>::template accessor<ro, ro> rk_dt_a) noexcept {

  auto [dt_mass, dt_energy, dt_mom] = RK<D>::mdcolex(m, rk_dt_a);

  double max_dm = 0.0;
  double max_du = 0.0;
  double max_de = 0.0;

  for(auto i : m.template cells<ax::x, dm::quantities>()) {
    max_dm = std::max(max_dm, std::abs(dt_mass(i)));
    max_du = std::max(max_du, std::abs(dt_mom(i)[0]));
    max_de = std::max(max_de, std::abs(dt_energy(i)));

    if(std::abs(dt_mass(i)) > 1e-12 || std::abs(dt_mom(i)[0]) > 1e-12 ||
       std::abs(dt_energy(i)) > 1e-12) {
      std::cout << "[RHS nonzero] i=" << i << " dt_mass=" << dt_mass(i)
                << " dt_mom_x=" << dt_mom(i)[0] << " dt_energy=" << dt_energy(i)
                << "\n";
    }
  }

  std::cout << "[RHS max] mass=" << max_dm << " mom_x=" << max_du
            << " energy=" << max_de << "\n";
}

// -----------------------------------------------------------------------------
// Geometry / validation helpers (device-callable)
// -----------------------------------------------------------------------------

template<std::size_t Dim>
FLECSI_INLINE_TARGET vec<Dim>
normalize_interface_normal(const vec<Dim> & n_left_face,
  const vec<Dim> & n_right_face,
  std::size_t face_axis) noexcept {

  vec<Dim> n_if = 0.5 * (n_left_face + n_right_face);
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

template<std::size_t Dim>
FLECSI_INLINE_TARGET bool
valid_face_state(double rho,
  const vec<Dim> & u,
  double p,
  double c,
  const vec<Dim> & ru,
  double re) noexcept {
  return std::isfinite(rho) && rho > 0.0 && std::isfinite(u.norm()) &&
         std::isfinite(p) && p > 0.0 && std::isfinite(c) && c > 0.0 &&
         std::isfinite(ru.norm()) && std::isfinite(re) && re > 0.0;
}

template<std::size_t Dim>
FLECSI_INLINE_TARGET flux_result<Dim>
physical_flux_from_state(std::size_t face_axis,
  const vec<Dim> & u,
  double p,
  const vec<Dim> & ru,
  double re) noexcept {

  flux_result<Dim> f;

  f.r_flux = ru[face_axis];

  f.ru_flux = ru[face_axis] * u;
  f.ru_flux[face_axis] += p;

  f.re_flux = (re + p) * u[face_axis];

  return f;
}

// -----------------------------------------------------------------------------
// HLLC flux wrapper (device-callable)
// -----------------------------------------------------------------------------

template<std::size_t Dim>
FLECSI_INLINE_TARGET flux_result<Dim>
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
  numerical_algorithms::wave_speed_est wse) noexcept {

  flux_result<Dim> result;

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
      pR, cR, f_r_H, numerical_algorithms::field_tag::rho, wse);

  result.ru_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, vec<Dim>>(face_axis,
      ruL, rhoL, uL, reL, pL, cL, f_ru_T, ruR, rhoR, uR, reR,
      pR, cR, f_ru_H, numerical_algorithms::field_tag::ru, wse);

  result.re_flux =
    numerical_algorithms::compute_HLLC_fluxes<Dim, double>(face_axis,
      reL, rhoL, uL, reL, pL, cL, f_rE_T, reR, rhoR, uR, reR,
      pR, cR, f_rE_H, numerical_algorithms::field_tag::re, wse);
  // clang-format on

  return result;
}

// This is the basic original contact GFM:
// pressure and velocity are copied from the real fluid at the ghost location.
// Set this true only if you intentionally want the tangential-slip extension.
inline constexpr bool gfm_preserve_tangential_jump = false;

template<std::size_t Dim>
FLECSI_INLINE_TARGET vec<Dim>
make_gfm_ghost_velocity(const vec<Dim> & n_if,
  const vec<Dim> & u_owner_extrapolated,
  const vec<Dim> & u_real_at_ghost_location) noexcept {

  if constexpr(Dim == 1) {
    return u_real_at_ghost_location;
  }
  else {
    if constexpr(gfm_preserve_tangential_jump) {
      const double un_owner = u_owner_extrapolated.dot(n_if);
      const double un_real = u_real_at_ghost_location.dot(n_if);

      // Owner tangential velocity + real normal velocity.
      return u_owner_extrapolated - un_owner * n_if + un_real * n_if;
    }
    else {
      // Basic GFM contact condition.
      return u_real_at_ghost_location;
    }
  }
}

// -----------------------------------------------------------------------------
// EOS-consistent helpers (templated on the device EOS dispatch type)
// -----------------------------------------------------------------------------

// EOS-consistent ghost state completion.
template<std::size_t Dim, typename EOS>
FLECSI_INLINE_TARGET void
compute_state_from_rho_u_p(const EOS & eos_disp,
  int mat_id,
  ghost_state<Dim> & ghost) noexcept {

  const auto mat = static_cast<std::size_t>(mat_id);
  const auto & eos = eos_disp[mat];

  ghost.e = hard::tasks::util::find_sie(eos, ghost.rho, ghost.p);
  ghost.c = eos.cRhoSie(ghost.rho, ghost.e);

  ghost.ru = ghost.rho * ghost.u;
  ghost.re = ghost.rho * ghost.e + 0.5 * ghost.rho * ghost.u.norm_squared();
}

template<std::size_t Dim, typename EOS>
FLECSI_INLINE_TARGET void
compute_state_from_rho_ru_re(const EOS & eos_disp,
  int mat_id,
  ghost_state<Dim> & ghost) noexcept {

  const auto mat = static_cast<std::size_t>(mat_id);
  const auto & eos = eos_disp[mat];

  ghost.u = ghost.ru / ghost.rho;
  ghost.e = ghost.re / ghost.rho - 0.5 * ghost.u.norm_squared();
  ghost.c = eos.cRhoSie(ghost.rho, ghost.e);
  ghost.p = eos.pRhoSie(ghost.rho, ghost.e);
}

template<std::size_t Dim>
FLECSI_INLINE_TARGET bool
valid_extp_q_for_gfm(const hard::tasks::util::extp_q<Dim> & q) noexcept {
  return std::isfinite(q.rho()) && q.rho() > 0.0 &&
         std::isfinite(q.u().norm()) && std::isfinite(q.p()) && q.p() > 0.0 &&
         std::isfinite(q.e()) && q.e() > 0.0;
}

template<std::size_t Dim>
FLECSI_INLINE_TARGET bool
valid_owner_extp_for_gfm(const hard::tasks::util::extp_q<Dim> & q) noexcept {
  return valid_extp_q_for_gfm<Dim>(q);
}

template<typename EOS>
FLECSI_INLINE_TARGET bool
valid_material_id(int mat_id, const EOS & eos_disp) noexcept {
  return mat_id >= 0 && static_cast<std::size_t>(mat_id) < eos_disp.size();
}

FLECSI_INLINE_TARGET bool
rel_close(double a,
  double b,
  double rtol = 1.0e-10,
  double atol = 1.0e-12) noexcept {
  return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

template<std::size_t Dim>
FLECSI_INLINE_TARGET double
sie_from_cons(double rho, const vec<Dim> & u, double re) noexcept {
  return (re - 0.5 * rho * u.norm_squared()) / rho;
}

template<std::size_t Dim, typename EOS>
FLECSI_INLINE_TARGET bool
eos_pair_equivalent_at_face(const EOS & eos_disp,
  int mL,
  int mR,
  double rhoL,
  const vec<Dim> & uL,
  double reL,
  double rhoR,
  const vec<Dim> & uR,
  double reR) noexcept {

  if(mL == mR) {
    return true;
  }

  if(!valid_material_id(mL, eos_disp) || !valid_material_id(mR, eos_disp)) {
    return false;
  }

  const double eL = sie_from_cons<Dim>(rhoL, uL, reL);
  const double eR = sie_from_cons<Dim>(rhoR, uR, reR);

  if(!std::isfinite(eL) || eL <= 0.0 || !std::isfinite(eR) || eR <= 0.0) {
    return false;
  }

  const auto & eosL = eos_disp[static_cast<std::size_t>(mL)];
  const auto & eosR = eos_disp[static_cast<std::size_t>(mR)];

  const double pLL = eosL.pRhoSie(rhoL, eL);
  const double pRL = eosR.pRhoSie(rhoL, eL);
  const double cLL = eosL.cRhoSie(rhoL, eL);
  const double cRL = eosR.cRhoSie(rhoL, eL);

  const double pLR = eosL.pRhoSie(rhoR, eR);
  const double pRR = eosR.pRhoSie(rhoR, eR);
  const double cLR = eosL.cRhoSie(rhoR, eR);
  const double cRR = eosR.cRhoSie(rhoR, eR);

  return rel_close(pLL, pRL) && rel_close(cLL, cRL) && rel_close(pLR, pRR) &&
         rel_close(cLR, cRR);
}

// Owner extrapolated-state sanity checks.
// Work-around to enable assert in device functions
template<std::size_t Dim>
FLECSI_INLINE_TARGET void
assert_owner_extp_state([[maybe_unused]] const char * label,
  [[maybe_unused]] const hard::tasks::util::extp_q<Dim> & q) noexcept {
#ifndef NDEBUG
  if(!(std::isfinite(q.rho()) && q.rho() > 0.0)) {
    std::printf("%s: invalid owner extrapolated rho\n", label);
    assert(false && "assert_owner_extp_state: invalid owner extrapolated rho");
  }

  if(!std::isfinite(q.u().norm())) {
    std::printf("%s: invalid owner extrapolated velocity\n", label);
    assert(
      false && "assert_owner_extp_state: invalid owner extrapolated velocity");
  }
#endif
}

// -----------------------------------------------------------------------------
// original GFM ghost construction
// -----------------------------------------------------------------------------
//
// qOwnerExt supplies the owner-material one-sided thermodynamic state.
// With the current interface, the available thermodynamic variable is rho.
// The ghost pressure and velocity are copied from the real fluid occupying the
// ghost location, then the state is completed with owner_mat EOS.
//

template<std::size_t Dim, typename EOS>
FLECSI_INLINE_TARGET ghost_state<Dim>
compute_original_gfm_ghost_state(const EOS & eos_disp,
  int owner_mat,
  std::size_t /* face_axis */,
  const vec<Dim> & n_if,
  const hard::tasks::util::extp_q<Dim> & qOwnerExt,
  const vec<Dim> & u_real_at_ghost_location,
  double p_real_at_ghost_location,
  const char * label) noexcept {

  assert_owner_extp_state<Dim>(label, qOwnerExt);

  assert(std::isfinite(p_real_at_ghost_location) &&
         p_real_at_ghost_location > 0.0 &&
         "compute_original_gfm_ghost_state: invalid copied real pressure");
  assert(std::isfinite(u_real_at_ghost_location.norm()) &&
         "compute_original_gfm_ghost_state: invalid copied real velocity");

  ghost_state<Dim> ghost;

  ghost.rho = qOwnerExt.rho();
  ghost.p = p_real_at_ghost_location;
  ghost.u =
    make_gfm_ghost_velocity<Dim>(n_if, qOwnerExt.u(), u_real_at_ghost_location);

  compute_state_from_rho_u_p<Dim>(eos_disp, owner_mat, ghost);

  return ghost;
}

// -----------------------------------------------------------------------------
// Riemann GFM ghost construction w/ entropy-minimum-preserving EoS
// -----------------------------------------------------------------------------
//
//
template<std::size_t Dim, typename EOS>
FLECSI_INLINE_TARGET ghost_state<Dim>
compute_riemann_gfm_ghost_state(const EOS & eos_disp,
  int owner_mat,
  int mL,
  int mR,
  const std::size_t face_axis,
  const vec<Dim> & n_if,
  const hard::tasks::util::extp_q<Dim> & qOwnerExt,
  double rho_real, // Real states are computed at ghost node location
  const vec<Dim> & u_real,
  double re_real,
  double p_real,
  double c_real,
  const char * label) noexcept {

  assert_owner_extp_state<Dim>(label, qOwnerExt);

  ghost_state<Dim> ghost;

  numerical_algorithms::wave_speed_est wse = numerical_algorithms::ClayTovWSE;
  ghost.rho = qOwnerExt.rho();
  ghost.p = qOwnerExt.p();
  ghost.u = qOwnerExt.u();
  compute_state_from_rho_u_p<Dim>(eos_disp, owner_mat, ghost);
  ghost_state<Dim> star_ghost;
  star_ghost.rho = numerical_algorithms::compute_HLLC_states<Dim, double>(n_if,
    ghost.rho,
    ghost.u,
    ghost.re,
    ghost.p,
    ghost.c, // Left state
    rho_real,
    u_real,
    re_real,
    p_real,
    c_real, // Right state
    numerical_algorithms::field_tag::rho,
    wse)
                     .first;

  star_ghost.ru = numerical_algorithms::compute_HLLC_states<Dim, vec<Dim>>(n_if,
    ghost.rho,
    ghost.u,
    ghost.re,
    ghost.p,
    ghost.c, // Left state
    rho_real,
    u_real,
    re_real,
    p_real,
    c_real, // Right state
    numerical_algorithms::field_tag::ru,
    wse)
                    .first;

  star_ghost.re = numerical_algorithms::compute_HLLC_states<Dim, double>(n_if,
    ghost.rho,
    ghost.u,
    ghost.re,
    ghost.p,
    ghost.c, // Left state
    rho_real,
    u_real,
    re_real,
    p_real,
    c_real, // Right state
    numerical_algorithms::field_tag::re,
    wse)
                    .first;

  // clang-format on

  compute_state_from_rho_ru_re<Dim>(eos_disp, owner_mat, star_ghost);

  return star_ghost;
}

// -----------------------------------------------------------------------------
// Process one face (device-callable)
// -----------------------------------------------------------------------------

template<std::size_t Dim, typename EOS>
FLECSI_INLINE_TARGET void
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
  const vec<Dim> & n_left_face,
  const vec<Dim> & n_right_face,
  const hard::tasks::util::extp_q<Dim> & qGhostR_L,
  const hard::tasks::util::extp_q<Dim> & qGhostL_R,
  const EOS & eos_disp,
  double & r_flux_out,
  vec<Dim> & ru_flux_out,
  double & re_flux_out,
  flux_result<Dim> & fluxes_L_out,
  flux_result<Dim> & fluxes_R_out,
  bool & is_cut_face) noexcept {

  const numerical_algorithms::wave_speed_est wse =
    numerical_algorithms::ClayTovWSE;

  is_cut_face = (mL != mR);

  if(!is_cut_face) {
    const auto fluxes = compute_all_HLLC_fluxes<Dim>(
      face_axis, rhoL, uL, pL, cL, ruL, reL, rhoR, uR, pR, cR, ruR, reR, wse);

    fluxes_L_out = fluxes;
    fluxes_R_out = fluxes;

    r_flux_out = fluxes.r_flux;
    ru_flux_out = fluxes.ru_flux;
    re_flux_out = fluxes.re_flux;

    return;
  }

  ghost_state<Dim> ghostR_L, ghostL_R;

  if constexpr(use_original_gfm) {

    const vec<Dim> n_if =
      normalize_interface_normal(n_left_face, n_right_face, face_axis);

    // Left-owner one-phase problem:
    //   real left material mL
    //   ghost right state with material mL EOS
    //   ghost p/u copied from real right state
    ghostR_L = compute_original_gfm_ghost_state<Dim>(
      eos_disp, mL, face_axis, n_if, qGhostR_L, uR, pR, "ghostR_L");

    // Right-owner one-phase problem:
    //   ghost left state with material mR EOS
    //   real right material mR
    //   ghost p/u copied from real left state
    ghostL_R = compute_original_gfm_ghost_state<Dim>(
      eos_disp, mR, face_axis, n_if, qGhostL_R, uL, pL, "ghostL_R");
  }
  else {
    ghostR_L = compute_riemann_gfm_ghost_state<Dim>(eos_disp,
      mL,
      mL,
      mR,
      face_axis,
      n_right_face,
      qGhostR_L,
      rhoR,
      uR,
      reR,
      pR,
      cR,
      "ghostR_L");
    ghostL_R = compute_riemann_gfm_ghost_state<Dim>(eos_disp,
      mR,
      mR,
      mL,
      face_axis,
      -1.0 * n_left_face,

      qGhostL_R,
      rhoL,
      uL,
      reL,
      pL,
      cL,
      "ghostL_R");
  }

  // clang-format off
  fluxes_L_out = compute_all_HLLC_fluxes<Dim>(face_axis,
    rhoL,
    uL,
    pL,
    cL,
    ruL,
    reL,
    ghostR_L.rho,
    ghostR_L.u,
    ghostR_L.p,
    ghostR_L.c,
    ghostR_L.ru,
    ghostR_L.re,
    wse);

  fluxes_R_out = compute_all_HLLC_fluxes<Dim>(face_axis,
    ghostL_R.rho,
    ghostL_R.u,
    ghostL_R.p,
    ghostL_R.c,
    ghostL_R.ru,
    ghostL_R.re,
    rhoR,
    uR,
    pR,
    cR,
    ruR,
    reR,
    wse);
  // clang-format on

  // Store left-owner flux in the face flux array.  The right/plus cell gets a
  // correction later so that it receives the right-owner flux.
  r_flux_out = fluxes_L_out.r_flux;
  ru_flux_out = fluxes_L_out.ru_flux;
  re_flux_out = fluxes_L_out.re_flux;
}

// -----------------------------------------------------------------------------
// Main task
// -----------------------------------------------------------------------------

struct compute_interface_fluxes_gfm_task : annotation_task_region {
  inline static const std::string name{"compute_interface_fluxes_gfm"};
};

template<std::size_t Dim>
void
compute_interface_fluxes_gfm(flecsi::exec::accelerator s,
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

  using common::tasks::utils::get_mdiota_policy;

  eos::eos_dispatch<flecsi::exec::accelerator> eos_disp(eos_mm);

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

  constexpr bool enable_gfm = true;

  // ---------------------------------------------------------------------------
  // 1D
  // ---------------------------------------------------------------------------
  if constexpr(Dim == 1) {
    flog_assert(
      face_axis == 0, "1D executable can only compute fluxes along X axis");

    const auto inv_dx = one_over_dx_i[0];
    const auto dx = 1.0 / inv_dx;

    // Face-flux pass (+ GFM cut correction to the plus/right cell).
    s.executor().forall(i, (m.template cells<ax::x, dm::corrector>())) {
      const double g_comp = (*g_acc)[0];

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

      const bool valid_material_pair = valid_material_id(mL_raw, eos_disp) &&
                                       valid_material_id(mR_raw, eos_disp);

      const bool face_state_valid_for_gfm =
        std::isfinite(rhoL) && rhoL > 0.0 && std::isfinite(rhoR) &&
        rhoR > 0.0 && std::isfinite(pL) && pL > 0.0 && std::isfinite(pR) &&
        pR > 0.0 && std::isfinite(cL) && cL > 0.0 && std::isfinite(cR) &&
        cR > 0.0 && std::isfinite(reL) && reL > 0.0 && std::isfinite(reR) &&
        reR > 0.0;

      const bool left_state_valid =
        valid_material_id(mL_raw, eos_disp) &&
        valid_face_state<Dim>(rhoL, u_left(i - 1), pL, cL, ru_left(i - 1), reL);

      const bool right_state_valid =
        valid_material_id(mR_raw, eos_disp) &&
        valid_face_state<Dim>(rhoR, u_right(i), pR, cR, ru_right(i), reR);

      // Physical boundary / invalid ghost-side guard.
      // For transmissive/outflow boundary behavior, use the valid interior
      // side.
      if(!left_state_valid || !right_state_valid) {
        flux_result<Dim> fb{};

        if(!left_state_valid && right_state_valid) {
          // Left physical boundary: use right/interior state.
          fb =
            physical_flux_from_state<Dim>(0, u_right(i), pR, ru_right(i), reR);
        }
        else if(left_state_valid && !right_state_valid) {
          // Right physical boundary: use left/interior state.
          fb = physical_flux_from_state<Dim>(
            0, u_left(i - 1), pL, ru_left(i - 1), reL);
        }
        else {
          // Both sides invalid: no meaningful flux available.
          assert(false && "Both sides of a face are invalid");
        }

        r_f(i) = fb.r_flux;
        ru_f(i) = fb.ru_flux;
        re_f(i) = fb.re_flux;

        return;
      }

      const bool qext_valid_for_cut =
        valid_owner_extp_for_gfm<Dim>(q_ext_plus(i)) &&
        valid_owner_extp_for_gfm<Dim>(q_ext_minus(i - 1));

      const bool equivalent_eos_pair =
        valid_material_pair && eos_pair_equivalent_at_face<Dim>(eos_disp,
                                 mL_raw,
                                 mR_raw,
                                 rhoL,
                                 u_left(i - 1),
                                 reL,
                                 rhoR,
                                 u_right(i),
                                 reR);

      const bool use_gfm_cut = enable_gfm && valid_material_pair &&
                               face_state_valid_for_gfm && (mL_raw != mR_raw) &&
                               qext_valid_for_cut && !equivalent_eos_pair;

      const int mL_eff = mL_raw;
      const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      flux_result<Dim> fluxes_L, fluxes_R;
      bool is_cut = false;

      // clang-format off
      process_single_face<Dim>(0,
        mL_eff,
        mR_eff,
        rhoL,
        u_left(i - 1),
        pL,
        cL,
        ru_left(i - 1),
        reL,
        rhoR,
        u_right(i),
        pR,
        cR,
        ru_right(i),
        reR,
        n_right(i - 1),     // left cell plus/right face
        n_left(i),          // right cell minus/left face
        q_ext_plus(i),
        q_ext_minus(i - 1),
        eos_disp,
        r_f(i),
        ru_f(i),
        re_f(i),
        fluxes_L,
        fluxes_R,
        is_cut);
      // clang-format on

      if(use_gfm_cut) {
        dt_mass_density(i) += inv_dx * (fluxes_R.r_flux - fluxes_L.r_flux);
        dt_momentum_density(i) +=
          inv_dx * (fluxes_R.ru_flux - fluxes_L.ru_flux);
        dt_total_energy_density(i) +=
          inv_dx * (fluxes_R.re_flux - fluxes_L.re_flux);
      }
    };

    // Flux-divergence pass.
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_mass_density(i) += inv_dx * (r_f(i) - r_f(i + 1));
      dt_momentum_density(i) += inv_dx * (ru_f(i) - ru_f(i + 1));
      dt_total_energy_density(i) += inv_dx * (re_f(i) - re_f(i + 1));
    };
  }

  // ---------------------------------------------------------------------------
  // 2D
  // ---------------------------------------------------------------------------
  else if constexpr(Dim == 2) {
    flog_assert(
      face_axis <= 1, "2D executable cannot compute fluxes along Z axis");

    const std::array<std::array<int, 2>, 2> left_offset = {{{-1, 0}, {0, -1}}};

    const auto offset = left_offset[face_axis];
    const auto inv_dx = one_over_dx_i[face_axis];

    auto mdpolicy = (face_axis == 0)
                      ? get_mdiota_policy(r_f,
                          m.template cells<ax::y, dm::quantities>(),
                          m.template cells<ax::x, dm::corrector>())
                      : get_mdiota_policy(r_f,
                          m.template cells<ax::y, dm::corrector>(),
                          m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy) {
      auto [j, i] = ji;

      const double g_comp = (*g_acc)[face_axis];

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

      const bool valid_material_pair = valid_material_id(mL_raw, eos_disp) &&
                                       valid_material_id(mR_raw, eos_disp);

      const bool face_state_valid_for_gfm =
        std::isfinite(rhoL) && rhoL > 0.0 && std::isfinite(rhoR) &&
        rhoR > 0.0 && std::isfinite(pL) && pL > 0.0 && std::isfinite(pR) &&
        pR > 0.0 && std::isfinite(cL) && cL > 0.0 && std::isfinite(cR) &&
        cR > 0.0 && std::isfinite(reL) && reL > 0.0 && std::isfinite(reR) &&
        reR > 0.0;

      const bool qext_valid_for_cut =
        valid_owner_extp_for_gfm<Dim>(q_ext_plus(i, j)) &&
        valid_owner_extp_for_gfm<Dim>(q_ext_minus(iL, jL));

      const bool use_gfm_cut = enable_gfm && valid_material_pair &&
                               face_state_valid_for_gfm && (mL_raw != mR_raw) &&
                               qext_valid_for_cut;

      const int mL_eff = mL_raw;
      const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      flux_result<Dim> fluxes_L, fluxes_R;
      bool is_cut = false;

      // clang-format off
      process_single_face<Dim>(face_axis,
        mL_eff,
        mR_eff,
        rhoL,
        u_left(iL, jL),
        pL,
        cL,
        ru_left(iL, jL),
        reL,
        rhoR,
        u_right(i, j),
        pR,
        cR,
        ru_right(i, j),
        reR,
        n_right(iL, jL),     // left cell plus face
        n_left(i, j),        // right cell minus face
        q_ext_plus(i, j),
        q_ext_minus(iL, jL),
        eos_disp,
        r_f(i, j),
        ru_f(i, j),
        re_f(i, j),
        fluxes_L,
        fluxes_R,
        is_cut);
      // clang-format on

      if(use_gfm_cut) {
        dt_mass_density(i, j) += inv_dx * (fluxes_R.r_flux - fluxes_L.r_flux);
        dt_momentum_density(i, j) +=
          inv_dx * (fluxes_R.ru_flux - fluxes_L.ru_flux);
        dt_total_energy_density(i, j) +=
          inv_dx * (fluxes_R.re_flux - fluxes_L.re_flux);
      }
    };

    auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
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
    };
  }

  // ---------------------------------------------------------------------------
  // 3D
  // ---------------------------------------------------------------------------
  else {
    const std::array<std::array<int, 3>, 3> left_offset = {
      {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}};

    const auto offset = left_offset[face_axis];
    const auto inv_dx = one_over_dx_i[face_axis];

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

    s.executor().forall(kji, mdpolicy) {
      auto [k, j, i] = kji;

      const double g_comp = (*g_acc)[face_axis];

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

      const bool valid_material_pair = valid_material_id(mL_raw, eos_disp) &&
                                       valid_material_id(mR_raw, eos_disp);

      const bool face_state_valid_for_gfm =
        std::isfinite(rhoL) && rhoL > 0.0 && std::isfinite(rhoR) &&
        rhoR > 0.0 && std::isfinite(pL) && pL > 0.0 && std::isfinite(pR) &&
        pR > 0.0 && std::isfinite(cL) && cL > 0.0 && std::isfinite(cR) &&
        cR > 0.0 && std::isfinite(reL) && reL > 0.0 && std::isfinite(reR) &&
        reR > 0.0;

      const bool qext_valid_for_cut =
        valid_owner_extp_for_gfm<Dim>(q_ext_plus(i, j, k)) &&
        valid_owner_extp_for_gfm<Dim>(q_ext_minus(iL, jL, kL));

      const bool use_gfm_cut = enable_gfm && valid_material_pair &&
                               face_state_valid_for_gfm && (mL_raw != mR_raw) &&
                               qext_valid_for_cut;

      const int mL_eff = mL_raw;
      const int mR_eff = use_gfm_cut ? mR_raw : mL_raw;

      flux_result<Dim> fluxes_L, fluxes_R;
      bool is_cut = false;

      // clang-format off
      process_single_face<Dim>(face_axis,
        mL_eff,
        mR_eff,
        rhoL,
        u_left(iL, jL, kL),
        pL,
        cL,
        ru_left(iL, jL, kL),
        reL,
        rhoR,
        u_right(i, j, k),
        pR,
        cR,
        ru_right(i, j, k),
        reR,
        n_right(iL, jL, kL),   // left cell plus face
        n_left(i, j, k),       // right cell minus face
        q_ext_plus(i, j, k),
        q_ext_minus(iL, jL, kL),
        eos_disp,
        r_f(i, j, k),
        ru_f(i, j, k),
        re_f(i, j, k),
        fluxes_L,
        fluxes_R,
        is_cut);
      // clang-format on

      if(use_gfm_cut) {
        dt_mass_density(i, j, k) +=
          inv_dx * (fluxes_R.r_flux - fluxes_L.r_flux);
        dt_momentum_density(i, j, k) +=
          inv_dx * (fluxes_R.ru_flux - fluxes_L.ru_flux);
        dt_total_energy_density(i, j, k) +=
          inv_dx * (fluxes_R.re_flux - fluxes_L.re_flux);
      }
    };

    auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
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
    };
  }
}

} // namespace hard::tasks::gfm

#endif // HARD_MODULE_GFM_INTERFACE_FLUXES_HH
