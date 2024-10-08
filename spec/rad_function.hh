#ifndef SPEC_RAD_FUNCTION_HH
#define SPEC_RAD_FUNCTION_HH

#include <cmath>
#include <flecsi/flog.hh>

namespace spec {

// Follwing structure contains definition and
// derivatives of a quartic polynomial
struct poly {
  // polynomial definition for root finding
  // see eq (37) of Moens et al. A&A 657, A81 (2022)
  // eq (37) has been multiplied by a1
  static double func(double x, double c2, double c1, double c0) {
    return c2 * x * x * x * x + c1 * x - c0;
  }

  // 1st derivative of the same polynomial
  static double funcprime(double x, double c2, double c1) {
    return 4.0 * c2 * x * x * x + c1;
  }

  // 2nd derivative of the same polynomial
  static double funcprime2(double x, double c2) {
    return 12.0 * c2 * x * x;
  }
};

struct rootFindingMethods {

  // Bisection Method
  struct Bisection {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // set tolerance
      double an = 0.0;
      double bn = bound;
      double cn;
      int iter = 0;
      const int max_iter = 3000; // safety maximum iteration limit

      while(iter < max_iter) {
        cn = 0.5 * (an + bn);
        double fa = poly::func(an, c2, c1, c0);
        double fc = poly::func(cn, c2, c1, c0);

        if(std::abs(fc) < tol) {
          return cn;
        }
        if(fa * fc < 0.0) {
          bn = cn;
        }
        else {
          an = cn;
        }
        iter++;
      }

      flog_fatal("Bisection Method failed!");

      return -1;
    };
  };

  // Newton-Raphson method
  struct NewtonRaphson {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // Tolerance
      double xn = bound * 0.5; // initial guess
      int iter = 0;
      const int max_iter = 3000; // safety maximum iteration limit

      while(iter < max_iter) {
        double fx = poly::func(xn, c2, c1, c0);
        double fprime = poly::funcprime(xn, c2, c1);

        if(std::abs(fx) < tol) {
          return xn;
        }

        xn = xn - fx / fprime;
        iter++;
      }

      flog_fatal("Newton-Raphson method failed!");
      return -1.0; // No convergence
    }; // get_root
  }; // NewtonRaphson

  // Secant method
  struct Secant {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // Tolerance
      double x0 = 0.0;
      double x1 = bound * 0.5;
      double fx0 = poly::func(x0, c2, c1, c0);
      double fx1 = poly::func(x1, c2, c1, c0);
      int iter = 0;
      const int max_iter = 3000; // Safety maximum iteration limit

      while(iter < max_iter) {
        if(std::abs(fx1) < tol) {
          return x1;
        }

        double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        x0 = x1;
        fx0 = fx1;
        x1 = x2;
        fx1 = poly::func(x1, c2, c1, c0);

        iter++;
      }

      flog_fatal("Secant method failed!");
      return -1.0; // No convergence
    }; // get_root
  }; // Secant

  // Halley's method
  struct HalleysMethod {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // Tolerance
      double xn = bound * 0.5; // Initial guess
      int iter = 0;
      const int max_iter = 3000; // Safety maximum iteration limit

      while(iter < max_iter) {
        double fx = poly::func(xn, c2, c1, c0);
        double fprime = poly::funcprime(xn, c2, c1);
        double fprime2 = poly::funcprime2(xn, c2);

        if(std::abs(fx) < tol) {
          return xn;
        }

        xn = xn - (2.0 * fx * fprime) / (2.0 * fprime * fprime - fx * fprime2);
        iter++;
      }

      flog_fatal("Halley's method failed!");
      return -1.0; // No convergence
    }; // get_root
  }; // HalleysMethod

  // Fixed point iteration
  struct FixedPointIteration {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // Tolerance
      double xn = bound * 0.5; // Initial guess
      int iter = 0;
      const int max_iter = 3000; // Safety maximum iteration limit

      while(iter < max_iter) {
        double gx = xn - poly::func(xn, c2, c1, c0); // Fixed-point function

        if(std::abs(gx - xn) < tol) {
          return xn;
        }

        xn = gx;
        iter++;
      }

      flog_fatal("Fixed point iteration method failed!");
      return -1.0; // No convergence
    }; // get_root
  }; // FixedPointIteration

  // Regula Falsi method (False Position)
  struct RegulaFalsi {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // Tolerance
      double an = 0.0;
      double bn = bound;
      double cn;
      int iter = 0;
      const int max_iter = 1000; // Safety maximum iteration limit

      while(iter < max_iter) {
        cn = an - poly::func(an, c2, c1, c0) * (bn - an) /
                    (poly::func(bn, c2, c1, c0) - poly::func(an, c2, c1, c0));

        if(std::abs(poly::func(cn, c2, c1, c1)) < tol) {
          return cn;
        }

        if(poly::func(an, c2, c1, c0) * poly::func(cn, c2, c1, c0) < 0.0) {
          bn = cn;
        }
        else {
          an = cn;
        }

        iter++;
      }

      flog_fatal("Regula Falsi method failed!");
      return -1.0; // No convergence
    }; // get_root
  }; // Regula Falsi Method

  // Newton-Raphson with bracketing method (Bounded Newton-Raphson)
  struct NewtonRaphsonBracketing {
    static double get_root(double c2, double c1, double c0, double bound) {
      const double tol = 1.0e-12; // tolerance
      double xn = bound * 0.5; // initial guess
      int iter = 0;
      const int max_iter = 1000; // safety maximum iteration limit

      while(iter < max_iter) {
        double fx = poly::func(xn, c2, c1, c0);
        double fprime = poly::funcprime(xn, c2, c1);

        if(std::abs(fx) < tol) {
          return xn;
        }

        double xnew = xn - fx / fprime;

        if(xnew < 0.0 || xnew > bound) {
          break;
        }

        xn = xnew;
        iter++;
      }

      // If not converged with Newton-Raphson, use Bisection as fallback
      double an = 0.0;
      double bn = bound;
      double cn;

      iter = 0;

      while(iter < max_iter) {
        cn = 0.5 * (an + bn);
        double fa = poly::func(an, c2, c1, c0);
        double fc = poly::func(cn, c2, c1, c0);

        if(std::abs(fc) < tol) {
          return cn;
        }

        if(fa * fc < 0.0) {
          bn = cn;
        }
        else {
          an = cn;
        }

        iter++;
      }

      flog_fatal("Newton-Raphson with bracketing method failed!");
      return -1.0; // No convergence
    }
  };

}; // rootfindingMethods
}; // namespace spec

#endif // SPEC_RAD_FUNCTION_HH
