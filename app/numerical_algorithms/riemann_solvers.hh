
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

  const double l_max = std::max(0.0, std::max(LmaxT, LmaxH));
  const double l_min = std::min(0.0, std::min(LminT, LminH));

  const double ll = l_max * l_min;

  return (l_max * flux_left - l_min * flux_right + ll * (right - left)) /
         (l_max - l_min);
}
} // namespace hard::numerical_algorithms
