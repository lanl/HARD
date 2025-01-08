#ifndef SPEC_LIMITER_HH
#define SPEC_LIMITER_HH

#include "utils.hh"
#include <cmath>

namespace spec::limiters {

struct genminmod {
  FLECSI_INLINE_TARGET static double limit(double qL, double qC, double qR) {
    const double GMMLTheta{1.5};
    const double sC{0.5 * (qR - qL)};
    const double sL{GMMLTheta * (qC - qL)};
    const double sR{GMMLTheta * (qR - qC)};
    return (sL * sR < 0.0
              ? 0.0
              : (sC > 0.0 ? 1.0 : -1.0) *
                  std::min(std::abs(sC), std::min(std::abs(sL), std::abs(sR))));
  } // operator()
}; // genminmod

struct minmod {
  FLECSI_INLINE_TARGET static std::array<double, 2>
  reconstruct(double u_minus, double u_0, double u_plus) {
    const double a = u_plus - u_0;
    const double b = u_0 - u_minus;
    const double slope = 0.5 * (std::copysign(1.0, a) + std::copysign(1.0, b)) *
                         std::min(std::abs(a), std::abs(b));
    return {{u_0 - 0.5 * slope, u_0 + 0.5 * slope}};
  }
};

struct ppm4 {

  FLECSI_INLINE_TARGET static std::array<double, 2> reconstruct(
    double u_minus_2,
    double u_minus,
    double u_0,
    double u_plus,
    double u_plus_2) {

    const double one_over_twelve{0.08333333333333333};

    // Unlimited 4-th order reconstruction
    double u_plus_12 =
      (7.0 * (u_0 + u_plus) - (u_plus_2 + u_minus)) * one_over_twelve;
    double u_minus_12 =
      (7.0 * (u_0 + u_minus) - (u_minus_2 + u_plus)) * one_over_twelve;

    // Limit the unlimited reconstruction values to lie between adjacent cell
    // averages.
    u_plus_12 = std::max(u_plus_12, std::min(u_0, u_plus));
    u_plus_12 = std::min(u_plus_12, std::max(u_0, u_plus));
    u_minus_12 = std::max(u_minus_12, std::min(u_0, u_minus));
    u_minus_12 = std::min(u_minus_12, std::max(u_0, u_minus));

    // Evaluate differences
    const double delta_u_plus = u_plus_12 - u_0;
    const double delta_u_minus = u_0 - u_minus_12;

    if(delta_u_plus * delta_u_minus < 0.0) {
      return {{u_0, u_0}};
    }
    else {

      if(std::abs(delta_u_plus) >= 2.0 * std::abs(delta_u_minus)) {
        u_plus_12 = u_0 + 2.0 * delta_u_minus;
      }
      if(std::abs(delta_u_minus) >= 2.0 * std::abs(delta_u_plus)) {
        u_minus_12 = u_0 - 2.0 * delta_u_plus;
      }

      return {{u_minus_12, u_plus_12}};
    }
  }
};

struct weno5z {
  FLECSI_INLINE_TARGET static std::array<double, 2> reconstruct(
    double u_minus_2,
    double u_minus,
    double u_0,
    double u_plus,
    double u_plus_2) {

    const double epsilon = 1.0e-15;

    using std::abs;
    using utils::sqr;

    const std::array<double, 3> beta{
      1.0833333333333333 * sqr(u_minus_2 - 2.0 * u_minus + u_0) +
        0.25 * sqr(u_minus_2 - 4.0 * u_minus + 3.0 * u_0),
      1.0833333333333333 * sqr(u_minus - 2.0 * u_0 + u_plus) +
        0.25 * sqr(u_plus - u_minus),
      1.0833333333333333 * sqr(u_plus_2 - 2.0 * u_plus + u_0) +
        0.25 * sqr(u_plus_2 - 4.0 * u_plus + 3.0 * u_0)};

    const double tau5{abs(beta[2] - beta[0])};

    const std::array<double, 3> epsilon_k{
      epsilon * (1.0 + abs(u_0) + abs(u_minus) + abs(u_minus_2)),
      epsilon * (1.0 + abs(u_0) + abs(u_minus) + abs(u_plus)),
      epsilon * (1.0 + abs(u_0) + abs(u_plus) + abs(u_plus_2))};

    const std::array<double, 3> nw_buffer{
      1.0 + sqr(tau5 / (beta[0] + epsilon_k[0])),
      1.0 + sqr(tau5 / (beta[1] + epsilon_k[1])),
      1.0 + sqr(tau5 / (beta[2] + epsilon_k[2]))};

    const std::array<double, 3> alpha_upper{
      nw_buffer[0], 10.0 * nw_buffer[1], 5.0 * nw_buffer[2]};
    const std::array<double, 3> alpha_lower{
      nw_buffer[2], 10.0 * nw_buffer[1], 5.0 * nw_buffer[0]};
    const double alpha_norm_upper =
      alpha_upper[0] + alpha_upper[1] + alpha_upper[2];
    const double alpha_norm_lower =
      alpha_lower[0] + alpha_lower[1] + alpha_lower[2];

    // reconstruction stencils
    const std::array<double, 3> recons_stencils_upper{
      0.375 * u_minus_2 - 1.25 * u_minus + 1.875 * u_0,
      -0.125 * u_minus + 0.75 * u_0 + 0.375 * u_plus,
      0.375 * u_0 + 0.75 * u_plus - 0.125 * u_plus_2};
    const std::array<double, 3> recons_stencils_lower{
      0.375 * u_plus_2 - 1.25 * u_plus + 1.875 * u_0,
      -0.125 * u_plus + 0.75 * u_0 + 0.375 * u_minus,
      0.375 * u_0 + 0.75 * u_minus - 0.125 * u_minus_2};

    // reconstructed solutions
    return {{(alpha_lower[0] * recons_stencils_lower[0] +
               alpha_lower[1] * recons_stencils_lower[1] +
               alpha_lower[2] * recons_stencils_lower[2]) /
               alpha_norm_lower,
      (alpha_upper[0] * recons_stencils_upper[0] +
        alpha_upper[1] * recons_stencils_upper[1] +
        alpha_upper[2] * recons_stencils_upper[2]) /
        alpha_norm_upper}};
  }
};

} // namespace spec::limiters

#endif // SPEC_LIMITER_HH
