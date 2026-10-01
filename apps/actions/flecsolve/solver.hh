#ifndef HARD_APPS_ACTIONS_FLECSOLVE_SOLVER_HH
#define HARD_APPS_ACTIONS_FLECSOLVE_SOLVER_HH

#include "flecsolve/solvers/solver_settings.hh"
#include "flecsolve/time-integrators/rk23.hh"
#include "flecsolve/util/scalar_ops.hh"
#include "spec/labels.hh"

namespace hard::time_integrator::solver {

namespace rk23 = flecsolve::time_integrator::rk23;
using rk23::parameters;

// Product operator between a scalar and a future, necessary
// for e.g., rk4 when doing 0.5 * dt.
using flecsolve::operator*;

// Fixed-step method.

// Create the vectors for the parameters
template<std::size_t nvecs>
flecsolve::work_factory<nvecs> make_work;

template<class O, class W>
using P = rk23::parameters<O, W, rk23::stepper_settings>;

struct advance_method : annotation_action_region {
  inline static const std::string name{"advance_method"};
};

template<spec::slv::solver S, class P>
struct stepper {
  stepper(P p) : params(std::move(p)) {}

  // Time integrator.
  template<class Current, class Output>
  void advance(flecsi::future<double> dt, Current & current, Output & output) {

    annotation_rguard<advance_method> guard;

    auto & rhs = params.get_operator();

    if constexpr(S == spec::slv::euler) {
      auto & [k1] = params.get_work();

      // k1 = F(U^n)
      rhs.apply(current, k1);

      // U^{n+1} = U^n + dt k1
      output.axpy(dt, k1, current);
    }
    else if constexpr(S == spec::slv::heun) {
      auto & [k1, k2, stage] = params.get_work();

      // k1 = F(U^n)
      rhs.apply(current, k1);

      // U_stage = U^n + dt k1
      stage.axpy(dt, k1, current);

      // k2 = F(U_stage)
      rhs.apply(stage, k2);

      // k1 <- 0.5 k1 + 0.5 k2
      k1.linear_sum(0.5, k1, 0.5, k2);

      // U^{n+1} = U^n + dt k1
      output.axpy(dt, k1, current);
    }
    else if constexpr(S == spec::slv::rk4) {
      auto & [k1, k2, k3, k4, next] = params.get_work();

      auto dt_1_2 = 0.5 * dt;
      auto dt_1_6 = (1. / 6.) * dt;

      // k1 = F(U^n)
      rhs.apply(current, k1);

      // k2 = F(U^n + dt/2 * k1)
      next.axpy(std::move(dt_1_2), k1, current);
      rhs.apply(next, k2);

      // k3 = F(U^n + dt/2 * k2)
      next.axpy(std::move(dt_1_2), k2, current);
      rhs.apply(next, k3);

      // k4 = F(U^n + dt * k3)
      next.axpy(dt, k3, current);
      rhs.apply(next, k4);

      // U^{n+1} = U^n + dt * next
      // With next = (k1 + 2 * k2 + 2 * k3 + k4) / 6
      next.linear_sum(1., k1, 2., k2);
      next.linear_sum(1., next, 2., k3);
      next.add(next, k4);
      output.axpy(std::move(dt_1_6), next, current);
    }
  }

protected:
  P params;
};

// Factory for the stepper so that the parameters can be deduced.
template<spec::slv::solver S, class P>
stepper<S, P>
make_stepper(P p) {
  return stepper<S, P>(std::move(p));
}

} // namespace hard::time_integrator::solver

#endif // HARD_APPS_ACTIONS_FLECSOLVE_SOLVER_HH
