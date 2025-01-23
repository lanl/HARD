
#pragma once

#include <cmath>
#include <cstddef>

#include <spec/types.hh>

namespace hard::numerical_algorithms {

//
// Returns HLL Riemann fluxes for the hydro part
//
template<std::size_t Dim>
FLECSI_INLINE_TARGET std::tuple<double, spec::vec<Dim>, double>
hll_hydro(
  // Evolved (conserved) variables
  const double mass_density_left,
  const double mass_density_right,
  const spec::vec<Dim> & momentum_density_left,
  const spec::vec<Dim> & momentum_density_right,
  const double total_energy_density_left,
  const double total_energy_density_right,
  // Fluxes
  const double flux_mass_density_left,
  const double flux_mass_density_right,
  const spec::vec<Dim> & flux_momentum_density_left,
  const spec::vec<Dim> & flux_momentum_density_right,
  const double flux_total_energy_density_left,
  const double flux_total_energy_density_right,
  // characteristic speeds from left and right states
  double min_char_speed_left,
  double max_char_speed_left,
  double min_char_speed_right,
  double max_char_speed_right) {

  std::tuple<double, spec::vec<Dim>, double> result{};

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
// Advection of radiation energy is perfomed using
// expressions 20 and 21 from
// Yang and Yuan, Publ. Astron. Soc. Japan 64, 69, 2012 August 25
FLECSI_INLINE_TARGET static double
advect_Erad(const double Erad_left,
  const double Erad_right,
  const double speed_left,
  const double speed_right) {

  const double speed_max =
    std::max(std::abs(speed_left), std::abs(speed_right));
  const double f_Erad_left{Erad_left * speed_left};
  const double f_Erad_right{Erad_right * speed_right};

  return 0.5 *
         (f_Erad_right + f_Erad_left - speed_max * (Erad_right - Erad_left));
}
} // namespace hard::numerical_algorithms
