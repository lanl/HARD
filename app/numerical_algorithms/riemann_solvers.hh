
#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include <spec/types.hh>

namespace flastro::numerical_algorithms {

//
// Returns HLL Riemann fluxes for the hydro part
//
template<std::size_t Dim>
FLECSI_INLINE_TARGET std::tuple<double, vec<Dim>, double>
hll_hydro(
  // Evolved (conserved) variables
  const double mass_density_left,
  const double mass_density_right,
  const vec<Dim> & momentum_density_left,
  const vec<Dim> & momentum_density_right,
  const double total_energy_density_left,
  const double total_energy_density_right,
  // Fluxes
  const double flux_mass_density_left,
  const double flux_mass_density_right,
  const vec<Dim> & flux_momentum_density_left,
  const vec<Dim> & flux_momentum_density_right,
  const double flux_total_energy_density_left,
  const double flux_total_energy_density_right,
  // characteristic speeds from left and right states
  double min_char_speed_left,
  double max_char_speed_left,
  double min_char_speed_right,
  double max_char_speed_right) {

  std::tuple<double, vec<Dim>, double> result{};

  const double lambda_max =
    std::max(0.0, std::max(max_char_speed_left, max_char_speed_right));
  const double lambda_min =
    std::min(0.0, std::min(min_char_speed_left, min_char_speed_right));

  const double lambdas_product = lambda_max * lambda_min;
  const double one_over_lambda_max_minus_lambda_min =
    1.0 / (lambda_max - lambda_min);

  get<0>(result) =
    (lambda_max * flux_mass_density_left -
      lambda_min * flux_mass_density_right +
      lambdas_product * (mass_density_right - mass_density_left)) *
    one_over_lambda_max_minus_lambda_min;

  // Note : vector operation is done here
  get<1>(result) =
    (lambda_max * flux_momentum_density_left -
      lambda_min * flux_momentum_density_right +
      lambdas_product * (momentum_density_right - momentum_density_left)) *
    one_over_lambda_max_minus_lambda_min;

  get<2>(result) = (lambda_max * flux_total_energy_density_left -
                     lambda_min * flux_total_energy_density_right +
                     lambdas_product * (total_energy_density_right -
                                         total_energy_density_left)) *
                   one_over_lambda_max_minus_lambda_min;

  return result;
}

//
// Returns HLLC Riemann fluxes for the hydro part
//
template<std::size_t Dim>
FLECSI_INLINE_TARGET std::tuple<double, vec<Dim>, double>
hllc_hydro(std::size_t face_axis,
  // Evolved (conserved) variables
  const double mass_density_left,
  const double mass_density_right,
  const vec<Dim> & momentum_density_left,
  const vec<Dim> & momentum_density_right,
  const double total_energy_density_left,
  const double total_energy_density_right,
  // Fluxes
  const double flux_mass_density_left,
  const double flux_mass_density_right,
  const vec<Dim> & flux_momentum_density_left,
  const vec<Dim> & flux_momentum_density_right,
  const double flux_total_energy_density_left,
  const double flux_total_energy_density_right,
  // Primitives
  const double pressure_left,
  const double pressure_right,
  const double velocity_left,
  const double velocity_right,
  // characteristic speeds from left and right states
  double min_char_speed_left,
  double max_char_speed_left,
  double min_char_speed_right,
  double max_char_speed_right) {

  std::tuple<double, vec<Dim>, double> result{};

  std::array<double, Dim> normal{0.0};
  normal.at(face_axis) = 1.0;

  const double lambda_max =
    std::max(0.0, std::max(max_char_speed_left, max_char_speed_right));
  const double lambda_min =
    std::min(0.0, std::min(min_char_speed_left, min_char_speed_right));

  const double lambdas_product = lambda_max * lambda_min;
  const double one_over_lambda_max_minus_lambda_min =
    1.0 / (lambda_max - lambda_min);

  // Compute lambda_star, the contact wave speed. (cf. Eq 10.37 of Toro2009)
  const double lambda_star =
    (pressure_right - pressure_left +
      momentum_density_left.get(face_axis) * (lambda_min - velocity_left) -
      momentum_density_right.get(face_axis) * (lambda_max - velocity_right)) /
    (mass_density_left * (lambda_min - velocity_left) -
      mass_density_right * (lambda_max - velocity_right));

  // Precompute common numerical factors for state variables in star region
  // (cf. Eq 10.39 of Toro2009)
  const double mult_left =
    (lambda_min - velocity_left) / (lambda_min - lambda_star);
  const double mult_right =
    (lambda_max - velocity_right) / (lambda_max - lambda_star);

  if(lambda_star >= 0.0) {
    get<0>(result) = flux_mass_density_left +
                     lambda_min * (mult_left - 1.0) * mass_density_left;
    for(size_t d = 0; d < Dim; ++d) {
      get<1>(result).get(d) =
        flux_momentum_density_left.get(d) +
        lambda_min * (mass_density_left * mult_left *
                         (lambda_star - velocity_left) * normal.at(d) +
                       (mult_left - 1.0) * momentum_density_left.get(d));
    }
    get<2>(result) =
      flux_total_energy_density_left +
      lambda_min *
        ((mult_left - 1.0) * (total_energy_density_left + pressure_left) +
          mult_left * mass_density_left * lambda_star *
            (lambda_star - velocity_left));
  }
  else {
    get<0>(result) = flux_mass_density_right +
                     lambda_max * (mult_right - 1.0) * mass_density_right;
    for(size_t d = 0; d < Dim; ++d) {
      get<1>(result).get(d) =
        flux_momentum_density_right.get(d) +
        lambda_max * (mass_density_right * mult_right *
                         (lambda_star - velocity_right) * normal.at(d) +
                       (mult_right - 1.0) * momentum_density_right.get(d));
    }
    get<2>(result) =
      flux_total_energy_density_right +
      lambda_max *
        ((mult_right - 1.0) * (total_energy_density_right + pressure_right) +
          mult_right * mass_density_right * lambda_star *
            (lambda_star - velocity_right));
  }

  return result;
}

// Advection of radiation energy is performed with the HLL flux.
// If both of the fluid velocity are very close to zero or they are almost same,
// then use LLF flux.
FLECSI_INLINE_TARGET static double
advect_Erad(const double Erad_left,
  const double Erad_right,
  const double v_left,
  const double v_right) {

  const double f_Erad_left{Erad_left * v_left};
  const double f_Erad_right{Erad_right * v_right};

  if(std::abs(v_left - v_right) < 1.0e-10) {
    const double abs_lambda_max =
      std::max(std::abs(speed_left), std::abs(speed_right));

    return 0.5 * (f_Erad_right + f_Erad_left -
                   abs_lambda_max * (Erad_right - Erad_left));
  }
  else {

    const double lambda_max = std::max(0.0, std::max(v_left, v_right));
    const double lambda_min = std::min(0.0, std::min(v_left, v_right));

    const double lambdas_product = lambda_max * lambda_min;
    const double one_over_lambda_max_minus_lambda_min =
      1.0 / (lambda_max - lambda_min);

    return (lambda_max * flux_total_energy_density_left -
             lambda_min * flux_total_energy_density_right +
             lambdas_product *
               (total_energy_density_right - total_energy_density_left)) *
           one_over_lambda_max_minus_lambda_min;
  }
}
} // namespace flastro::numerical_algorithms
