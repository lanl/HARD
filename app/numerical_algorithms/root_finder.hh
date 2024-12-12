
#pragma once
#include "../rad_function.hh"
#include <assert.h>

namespace flastro::numerical_algorithms {

namespace root_finder {

constexpr double tol = 1.0e-12;
constexpr size_t max_iter = 1000;

struct Bisection {
  static double
  get_root(double c4, double c1, double c0, double upper_bracket_bound);
};

struct NewtonRaphson {
  static double
  get_root(double c4, double c1, double c0, double upper_bracket_bound);
};

struct Secant {
  static double
  get_root(double c4, double c1, double c0, double upper_bracket_bound);
};

// struct Halleys {
//   FLECSI_INLINE_TARGET static double
//   get_root(double c4, double c1, double c0, double upper_bracket_bound);
// };

FLECSI_INLINE_TARGET double
halleys_get_root(double c4, double c1, double c0, double upper_bracket_bound) {
  using namespace flastro::energy_polynomial;

  double xn = upper_bracket_bound * 0.5; // Initial guess
  double dx;

  for(std::size_t i = 0; i < root_finder::max_iter; ++i) {
    double fx = func(xn, c4, c1, c0);
    double fprime = funcprime(xn, c4, c1);
    double fprime2 = funcprime2(xn, c4);

    dx = -(2.0 * fx * fprime) / (2.0 * fprime * fprime - fx * fprime2);

    if((std::abs(fx) < tol) or (std::abs(dx) < tol * xn)) {
      assert((xn >= 0.0) && "Halley's method failed"); // TODO: wrap in IFDEBUG
      return xn;
    }

    xn = xn + dx;
  }

  assert(false && "Halley's method failed to converge");
  return -1.0; // No convergence
}

struct FixedPointIteration {
  static double
  get_root(double c4, double c1, double c0, double upper_bracket_bound);
};

struct RegulaFalsi {
  static double
  get_root(double c4, double c1, double c0, double upper_bracket_bound);
};

} // namespace root_finder
} // namespace flastro::numerical_algorithms
