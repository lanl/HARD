#ifndef HARD_APPS_ACTIONS_FLECSOLVE_SOLVER_HH
#define HARD_APPS_ACTIONS_FLECSOLVE_SOLVER_HH

#include "flecsolve/solvers/solver_settings.hh"
#include "flecsolve/time-integrators/base.hh"
#include "flecsolve/time-integrators/parameters.hh"
#include "spec/labels.hh"

namespace hard::time_integrator::solver {

// Fixed-step method.

// This setting exists so the local integrator follows the same general
// interface as the FleCSolve time integrators.
struct settings : flecsolve::time_integrator::base_settings {};

// Combine the time-integrator settings, RHS operator, and work vectors.
template<class Op, class Work>
struct parameters : flecsolve::time_integrator::parameters<settings, Op, Work> {

  using base = flecsolve::time_integrator::parameters<settings, Op, Work>;

  template<class W>
  parameters(const settings & s, flecsolve::op::handle<Op> op, W && work)
    : base(s, op, std::forward<W>(work)) {}
};

// Class-template argument deduction for parameters.
template<class Op, class Work>
parameters(const settings &, flecsolve::op::handle<Op>, Work &&)
  -> parameters<Op, Work>;

// Create the vectors for the parameters
template<std::size_t nvecs>
flecsolve::work_factory<nvecs> make_work;

// Time integrator.
template<spec::slv::solver S, class Op, class Work>
struct integrator : flecsolve::time_integrator::base<parameters<Op, Work>> {

  using parameter_type = parameters<Op, Work>;
  using base = flecsolve::time_integrator::base<parameter_type>;

  using base::assert_can_advance;
  using base::current_dt;
  using base::params;

  explicit integrator(parameter_type p) : base(std::move(p)) {}

  template<class Current, class Output>
  void advance(double dt, Current & current, Output & output) {
    assert_can_advance();
    current_dt = dt;

    auto & rhs = params.get_operator();

    if constexpr(S == spec::slv::euler) {
      auto & [k1] = params.work;

      // k1 = F(U^n)
      rhs.apply(current, k1);

      // U^{n+1} = U^n + dt k1
      output.axpy(dt, k1, current);
    }
    else if constexpr(S == spec::slv::heun) {
      auto & [k1, k2, stage] = params.work;

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
      auto & [k1, k2, k3, k4, next] = params.work;

      // k1 = F(U^n)
      rhs.apply(current, k1);

      // k2 = F(U^n + dt/2 * k1)
      next.axpy(0.5 * dt, k1, current);
      rhs.apply(next, k2);

      // k3 = F(U^n + dt/2 * k2)
      next.axpy(0.5 * dt, k2, current);
      rhs.apply(next, k3);

      // k4 = F(U^n + dt * k3)
      next.axpy(dt, k3, current);
      rhs.apply(next, k4);

      // U^{n+1} = U^n + dt * next
      // With next = (k1 + 2 * k2 + 2 * k3 + k4) / 6
      next.linear_sum(1., k1, 2., k2);
      next.linear_sum(1., next, 2., k3);
      next.add(next, k4);
      output.axpy(dt / 6., next, current);
    }
  }
};

// Factory for the integrator so that Op and Work can be deduced.
template<spec::slv::solver S, class Op, class Work>
integrator<S, Op, Work>
make_integrator(parameters<Op, Work> p) {
  return integrator<S, Op, Work>(std::move(p));
}

} // namespace hard::time_integrator::solver

#endif // HARD_APPS_ACTIONS_FLECSOLVE_SOLVER_HH
