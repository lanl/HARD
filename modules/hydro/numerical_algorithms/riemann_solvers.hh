#ifndef HARD_MODULE_HYDRO_RIEMANN_SOLVERS_HH
#define HARD_MODULE_HYDRO_RIEMANN_SOLVERS_HH

#include <cmath>
#include <tuple>

#include <../modules/spec/types.hh>

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
compute_S_star(const double u_l,
  const double u_r,
  const double rho_l,
  const double rho_r,
  const double p_l,
  const double p_r,
  const double c_l,
  const double c_r) {

  const double sqrt_rho_l = std::sqrt(rho_l);
  const double sqrt_rho_r = std::sqrt(rho_r);
  const double denom = sqrt_rho_l + sqrt_rho_r;

  const double u_hat = (u_l * sqrt_rho_l + u_r * sqrt_rho_r) / denom;
  const double c_hat =
    std::sqrt((c_l * c_l * sqrt_rho_l + c_r * c_r * sqrt_rho_r) / denom +
              (sqrt_rho_l * sqrt_rho_r) / (2.0 * denom * denom) * (u_r - u_l) *
                (u_r - u_l));

  const double s_l = std::min(u_l - c_l, u_hat - c_hat);
  const double s_r = std::max(u_r + c_r, u_hat + c_hat);

  const double s_star =
    (p_r - p_l + rho_l * u_l * (s_l - u_l) - rho_r * u_r * (s_r - u_r)) /
    (rho_l * (s_l - u_l) - rho_r * (s_r - u_r));

  return std::make_tuple(s_l, s_r, s_star);
}

// Compute U* for HLLC flux computation
template<std::size_t Dim, typename T>
FLECSI_INLINE_TARGET std::pair<T, T>
compute_U_star(const double u_l_n,
  const double u_r_n,
  const vec<Dim> vel_l_star,
  const vec<Dim> vel_r_star,
  const double rho_l,
  const double rho_r,
  const double e_l,
  const double e_r,
  const double p_l,
  const double p_r,
  const double s_l,
  const double s_r,
  const double s_star,
  const std::string & var_name) {

  T u_l_star(0.0);
  T u_r_star(0.0);

  double const fac_l = rho_l * (s_l - u_l_n) / (s_l - s_star);
  double const fac_r = rho_r * (s_r - u_r_n) / (s_r - s_star);

  if constexpr(std::is_same_v<T, double>) {
    if(var_name == "rho") {
      u_l_star = fac_l;
      u_r_star = fac_r;
    }
    else if(var_name == "E") {
      u_l_star =
        fac_l * (e_l / rho_l +
                  (s_star - u_l_n) * (s_star + p_l / (rho_l * (s_l - u_l_n))));
      u_r_star =
        fac_r * (e_r / rho_r +
                  (s_star - u_r_n) * (s_star + p_r / (rho_r * (s_r - u_r_n))));
    }
    else {
      assert(false && "Invalid var_name for double star state");
    }
  }
  else if constexpr(std::is_same_v<T, vec<Dim>>) {
    if(var_name == "rhou") {
      u_l_star = fac_l * vel_l_star;
      u_r_star = fac_r * vel_r_star;
    }
    else {

      assert(false && "Invalid var_name for double star state");
    }
  }
  return {u_l_star, u_r_star};
}

// Computes F* for HLLC flux computation
template<typename T>
FLECSI_INLINE_TARGET T
compute_F_star(const T f_l,
  const T f_r,
  const T q_l,
  const T q_r,
  const T u_l_star,
  const T u_r_star,
  const double s_l,
  const double s_r,
  const double s_star) {

  // Star fluxes
  const T f_l_star = f_l + s_l * (u_l_star - q_l);
  const T f_r_star = f_r + s_r * (u_r_star - q_r);

  // Final flux selection
  if(0.0 <= s_l) {
    return f_l;
  }
  else if(s_l <= 0.0 && 0.0 <= s_star) {
    return f_l_star;
  }
  else if(s_star <= 0.0 && 0.0 <= s_r) {
    return f_r_star;
  }
  else {
    return f_r;
  }
}

// computes and applies  HLLC fluxes to a homogeneous conservation equation
template<std::size_t Dim, typename T>
FLECSI_INLINE_TARGET T
compute_HLLC_fluxes(std::size_t fa,
  const T q_l,
  const double rho_l,
  const vec<Dim> u_l,
  const double e_l,
  const double p_l,
  const double c_l,
  const T f_l,
  const T q_r,
  const double rho_r,
  const vec<Dim> u_r,
  const double e_r,
  const double p_r,
  const double c_r,
  const T f_r,
  const std::string & var_name) {

  double u_l_n = 0.0, u_r_n = 0.0;

  if(fa == 0) {
    u_l_n = u_l.x();
    u_r_n = u_r.x();
  }
  if constexpr(Dim > 1)
    if(fa == 1) {
      u_l_n = u_l.y();
      u_r_n = u_r.y();
    }
  if constexpr(Dim > 2)
    if(fa == 2) {
      u_l_n = u_l.z();
      u_r_n = u_r.z();
    }

  auto [s_l, s_r, s_star] =
    compute_S_star(u_l_n, u_r_n, rho_l, rho_r, p_l, p_r, c_l, c_r);

  vec<Dim> vel_l_star = u_l, vel_r_star = u_r;
  if(fa == 0) {
    vel_l_star.x() = s_star;
    vel_r_star.x() = s_star;
  }
  if constexpr(Dim > 1)
    if(fa == 1) {
      vel_l_star.y() = s_star;
      vel_r_star.y() = s_star;
    }
  if constexpr(Dim > 2)
    if(fa == 2) {
      vel_l_star.z() = s_star;
      vel_r_star.z() = s_star;
    }

  // clang-format off
  auto [u_l_star, u_r_star] = compute_U_star<Dim, T>(u_l_n, u_r_n,
    vel_l_star, vel_r_star, rho_l, rho_r, e_l, e_r, p_l, p_r, s_l, s_r, s_star, var_name);
  // clang-format on

  return compute_F_star<T>(
    f_l, f_r, q_l, q_r, u_l_star, u_r_star, s_l, s_r, s_star);
}
} // namespace hard::numerical_algorithms

#endif // HARD_MODULE_HYDRO_RIEMANN_SOLVERS_HH
