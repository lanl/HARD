#ifndef HARD_APPS_ACTIONS_FLECSOLVE_HEUN_HH
#define HARD_APPS_ACTIONS_FLECSOLVE_HEUN_HH

#include "flecsolve/solvers/solver_settings.hh"
#include "flecsolve/time-integrators/base.hh"
#include "flecsolve/time-integrators/parameters.hh"

#include <cstddef>
#include <utility>

namespace hard::time_integrator::heun {

/*
  Fixed-step Heun method.

  This setting exists so the local integrator follows the same general
  interface as the FleCSolve time integrators.
*/
struct settings : flecsolve::time_integrator::base_settings {
  bool use_fixed_dt = true;
};

/*
  Combine the time-integrator settings, RHS operator, and work vectors.
*/
template<class Op, class Work>
struct parameters : flecsolve::time_integrator::parameters<settings, Op, Work> {

  using base = flecsolve::time_integrator::parameters<settings, Op, Work>;

  template<class W>
  parameters(const settings & s, flecsolve::op::handle<Op> op, W && work)
    : base(s, op, std::forward<W>(work)) {}
};

/*
  Class-template argument deduction for parameters.
*/
template<class Op, class Work>
parameters(const settings &, flecsolve::op::handle<Op>, Work &&)
  -> parameters<Op, Work>;

/*
  Heun requires:

    k1    = first RHS evaluation
    k2    = second RHS evaluation
    stage = intermediate state
*/
enum workvecs : std::size_t { k1, k2, stage, nvecs };

/*
  Create the required work vectors from a prototype vector.
*/
static inline flecsolve::work_factory<workvecs::nvecs> make_work;

/*
  Optional topology-work alias, following the FleCSolve integrator pattern.
*/
template<std::size_t Version = 0>
using topo_work = flecsolve::topo_work_base<workvecs::nvecs, Version>;

/*
  Heun time integrator.

  Given

      dU/dt = F(U),

  perform

      k1        = F(U^n)
      U_stage   = U^n + dt k1
      k2        = F(U_stage)
      U^{n+1}   = U^n + dt/2 (k1 + k2)
*/
template<class Op, class Work>
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

    auto & [k1, k2, stage] = params.work;

    /*
      Stage 1:

          k1 = F(U^n)
    */
    rhs.apply(current, k1);

    /*
      Intermediate state:

          U_stage = U^n + dt k1
    */
    stage.axpy(dt, k1, current);

    /*
      Stage 2:

          k2 = F(U_stage)
    */
    rhs.apply(stage, k2);

    /*
      Average the two slopes:

          k1 <- 0.5 k1 + 0.5 k2
    */
    k1.linear_sum(0.5, k1, 0.5, k2);

    /*
      Final Heun solution:

          U^{n+1} = U^n + dt k1

      where k1 now contains the averaged slope.
    */
    output.axpy(dt, k1, current);
  }
};

/*
  Class-template argument deduction for integrator.
*/
template<class Op, class Work>
integrator(parameters<Op, Work>) -> integrator<Op, Work>;

} // namespace hard::time_integrator::heun

#endif // HARD_APPS_ACTIONS_FLECSOLVE_HEUN_HH