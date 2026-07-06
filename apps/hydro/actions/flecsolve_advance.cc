#include "state.hh"

#include "../../actions/flecsolve/hydro_flecsolve.hh"

#include <spec/runtime.hh>
#include <flecsi/flog.hh>

namespace hard {

template<std::size_t D>
struct hydro_flecsolve_advance {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    /*
      FleCSolve-style hydro RHS operator.
    */
    actions::flecsolve_hydro::rhs_params<D> params{&s, &sc};
    actions::flecsolve_hydro::rhs_operator<D> rhs(params);

    /*
      Read dt from HARD's timestep field.

      The old HARD update tasks used s.dt directly inside tasks.
      The FleCSolve vector algebra needs dt as a scalar.
    */
    const double dt =
      sc.template execute<actions::flecsolve_hydro::read_dt_task>(
          flecsi::exec::on, s.dt(*s.gt))
        .get();

    /*
      Vector views:

          U  -> current conserved state, s.cons.hydro
          Un -> saved beginning-of-step state, s.rk_n
          R1 -> first RHS, s.rk_dt1
          R2 -> second RHS, s.rk_dt2
    */
    auto U = actions::flecsolve_hydro::make_state_vector(s);
    auto Un = actions::flecsolve_hydro::make_saved_state_vector(s);
    auto R1 = actions::flecsolve_hydro::make_rhs_vector(s);
    auto R2 = actions::flecsolve_hydro::make_rhs2_vector(s);

    /*
      Stage 1:

          R1 = F(U^n)
    */
    rhs.apply(U, R1);

    /*
      Stage 1 update:

          U* = U^n + dt R1

      This replaces update_rk_stage_1.
    */
    U.axpy(dt, R1, Un);

    actions::flecsolve_hydro::recover_primitives_and_apply_boundaries(s, sc);

    /*
      Stage 2:

          R2 = F(U*)
    */
    rhs.apply(U, R2);

    /*
      Average the two RK slopes:

          R1 <- 0.5 R1 + 0.5 R2
    */
    R1.linear_sum(0.5, R1, 0.5, R2);

    /*
      Final Heun update:

          U^{n+1} = U^n + dt R1

      Since R1 now stores the averaged slope, this is equivalent to:

          U^{n+1} = U^n + 0.5 dt (F(U^n) + F(U*))
    */
    U.axpy(dt, R1, Un);

    actions::flecsolve_hydro::recover_primitives_and_apply_boundaries(s, sc);
  }
};

static const auto flecsolve_advance_action =
  spec::register_action<control, state, hydro_flecsolve_advance, cp::rk_stage_1>();

} // namespace hard