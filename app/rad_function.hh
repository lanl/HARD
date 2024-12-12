
#pragma once

#include <spec/utils.hh>

namespace flastro {

//
// Following struct contains the definition and derivatives of a quartic
// polynomial that needs to be solved for updating gas and radiation energy
// densities.
//
// We use Eq (37) of Moens et al. A&A 657, A81 (2022), but the whole polynomial
// multiplied by the coefficient a1, to avoid a problematic behavior when the
// opacity approaches close to zero.
//
// f(x) = a1 * x^4 + (1 + a2) * x - (1 + a2) * en - a2 * En
//     == c4 * x^4 + c1 * x + c0
//
namespace energy_polynomial {

FLECSI_INLINE_TARGET double
func(double x, double c4, double c1, double c0) {

  return c4 * spec::utils::sqr(spec::utils::sqr(x)) + c1 * x + c0;
}

// 1st derivative of the same polynomial
FLECSI_INLINE_TARGET double
funcprime(double x, double c4, double c1) {
  return 4.0 * c4 * x * spec::utils::sqr(x) + c1;
}

// 2nd derivative of the same polynomial
FLECSI_INLINE_TARGET double
funcprime2(double x, double c4) {
  return 12.0 * c4 * spec::utils::sqr(x);
}
}; // namespace energy_polynomial

} // namespace flastro
