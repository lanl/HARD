
#pragma once

#include <cmath>

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
} // namespace hard::numerical_algorithms
