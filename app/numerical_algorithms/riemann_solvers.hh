
#pragma once

#include <cmath>
#include <tuple>

#include <spec/types.hh>

namespace hard::numerical_algorithms {

// Applies hll fluxes to a homogeneous conservation equation
template<typename T>
FLECSI_INLINE_TARGET T
advect_conserved(const T left,
  const T right,
  const T flux_left,
  const T flux_right,
  const double LminT,
  const double LmaxT,
  const double LminH,
  const double LmaxH) {

  const double Lmax = std::max(0.0, std::max(LmaxT, LmaxH));
  const double Lmin = std::min(0.0, std::min(LminT, LminH));

  return (Lmax * flux_left - Lmin * flux_right + Lmax * Lmin * (right - left)) /
         (Lmax - Lmin);
}

// Compute S* for HLLC flux computation
FLECSI_INLINE_TARGET std::tuple<double, double, double>
compute_S_star(const double uL,
  const double uR,
  const double rhoL,
  const double rhoR,
  const double pL,
  const double pR,
  const double cL,
  const double cR) {

  const double sqrt_rhoL = std::sqrt(rhoL);
  const double sqrt_rhoR = std::sqrt(rhoR);
  const double denom = sqrt_rhoL + sqrt_rhoR;

  const double u_hat = (uL * sqrt_rhoL + uR * sqrt_rhoR) / denom;
  const double c_hat = std::sqrt(
    (cL * cL * sqrt_rhoL + cR * cR * sqrt_rhoR) / denom +
    (sqrt_rhoL * sqrt_rhoR) / (2.0 * denom * denom) * (uR - uL) * (uR - uL));

  const double SL = std::min(uL - cL, u_hat - c_hat);
  const double SR = std::max(uR + cR, u_hat + c_hat);

  const double S_star =
    (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
    (rhoL * (SL - uL) - rhoR * (SR - uR));

  return std::make_tuple(SL, SR, S_star);
}

// Compute U* for HLLC flux computation
template<std::size_t Dim, typename T>
FLECSI_INLINE_TARGET std::pair<T, T>
compute_U_star(const double uL_n,
  const double uR_n,
  const vec<Dim> ul_star,
  const vec<Dim> ur_star,
  const double rhoL,
  const double rhoR,
  const double EL,
  const double ER,
  const double pL,
  const double pR,
  const double SL,
  const double SR,
  const double S_star,
  const std::string & var_name) {

  T UL_star(0.0);
  T UR_star(0.0);

  double const facL = rhoL * (SL - uL_n) / (SL - S_star);
  double const facR = rhoR * (SR - uR_n) / (SR - S_star);

  if constexpr(std::is_same_v<T, double>) {
    if(var_name == "rho") {
      UL_star = facL;
      UR_star = facR;
    }
    else if(var_name == "E") {
      UL_star = facL * (EL / rhoL + (S_star - uL_n) *
                                      (S_star + pL / (rhoL * (SL - uL_n))));
      UR_star = facR * (ER / rhoR + (S_star - uR_n) *
                                      (S_star + pR / (rhoR * (SR - uR_n))));
    }
    else {
      assert(false && "Invalid var_name for double star state");
    }
  }
  else if constexpr(std::is_same_v<T, vec<Dim>>) {
    if(var_name == "rhou") {
      UL_star = facL * ul_star;
      UR_star = facR * ur_star;
    }
    else {

      assert(false && "Invalid var_name for double star state");
    }
  }
  return {UL_star, UR_star};
}

// Computes F* for HLLC flux computation
template<typename T>
FLECSI_INLINE_TARGET T
compute_F_star(const T FL,
  const T FR,
  const T QL,
  const T QR,
  const T UL_star,
  const T UR_star,
  const double SL,
  const double SR,
  const double S_star) {

  // Star fluxes
  const T FL_star = FL + SL * (UL_star - QL);
  const T FR_star = FR + SR * (UR_star - QR);

  // Final flux selection
  if(0.0 <= SL) {
    return FL;
  }
  else if(SL <= 0.0 && 0.0 <= S_star) {
    return FL_star;
  }
  else if(S_star <= 0.0 && 0.0 <= SR) {
    return FR_star;
  }
  else {
    return FR;
  }
}

// computes and applies  HLLC fluxes to a homogeneous conservation equation
template<std::size_t Dim, typename T>
FLECSI_INLINE_TARGET T
compute_HLLC_fluxes(std::size_t fa,
  const T QL,
  const double rhoL,
  const vec<Dim> uL,
  const double EL,
  const double pL,
  const double cL,
  const T FL,
  const T QR,
  const double rhoR,
  const vec<Dim> uR,
  const double ER,
  const double pR,
  const double cR,
  const T FR,
  const std::string & var_name) {

  double uL_n = 0.0, uR_n = 0.0;

  if(fa == 0) {
    uL_n = uL.x();
    uR_n = uR.x();
  }
  if constexpr(Dim > 1)
    if(fa == 1) {
      uL_n = uL.y();
      uR_n = uR.y();
    }
  if constexpr(Dim > 2)
    if(fa == 2) {
      uL_n = uL.z();
      uR_n = uR.z();
    }

  auto [SL, SR, S_star] =
    compute_S_star(uL_n, uR_n, rhoL, rhoR, pL, pR, cL, cR);

  vec<Dim> ul_star = uL, ur_star = uR;
  if(fa == 0) {
    ul_star.x() = S_star;
    ur_star.x() = S_star;
  }
  if constexpr(Dim > 1)
    if(fa == 1) {
      ul_star.y() = S_star;
      ur_star.y() = S_star;
    }
  if constexpr(Dim > 2)
    if(fa == 2) {
      ul_star.z() = S_star;
      ur_star.z() = S_star;
    }

  // clang-format off
  auto [UL_star, UR_star] = compute_U_star<Dim, T>(uL_n, uR_n,
    ul_star, ur_star, rhoL, rhoR, EL, ER, pL,pR, SL, SR, S_star, var_name);
  // clang-format on

  return compute_F_star<T>(FL, FR, QL, QR, UL_star, UR_star, SL, SR, S_star);
}
} // namespace hard::numerical_algorithms
