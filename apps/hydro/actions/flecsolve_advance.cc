#include "state.hh"

#include "../../actions/flecsolve/heun.hh"
#include "../../actions/flecsolve/hydro_flecsolve.hh"

#include <flecsi/flog.hh>
#include <spec/runtime.hh>

// #include "flecsolve/time-integrators/rk23.hh"

namespace hard {

namespace heun = time_integrator::heun;

template<std::size_t D>
struct hydro_flecsolve_advance {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    /*
    Read the time-step selected by HARD.
    */
    const double dt =
      sc.template execute<actions::flecsolve_hydro::read_dt_task>(
          flecsi::exec::on, s.dt(*s.gt))
        .get();

    if(dt <= 0.0) {
      return;
    }

    /*
      HARD hydro RHS wrapped as a FleCSolve operator.

      Do not wrap this in flecsolve::op::core. The core wrapper tries to
      subset anonymous FleCSolve variables, but our hydro multivector uses
      hard::actions::flecsolve_hydro::var.
    */
    actions::flecsolve_hydro::rhs_params<D> params{&s, &sc};
    actions::flecsolve_hydro::rhs_operator<D> F(params);

    /*
    U is a view of the current HARD conserved state.
    Uout uses the existing rk_n storage as the destination for the completed
    heun step.
    */
    auto U = actions::flecsolve_hydro::make_state_vector(s);
    auto Uout = actions::flecsolve_hydro::make_saved_state_vector(s);

    /*
    Configure one fixed Heun step over the timestep selected by HARD
    */
    heun::settings settings{};

    settings.initial_time = 0.0;
    settings.final_time = dt;
    settings.max_steps = 1;

    settings.initial_dt = dt;
    settings.min_dt = dt;
    settings.max_dt = dt;

    settings.use_fixed_dt = true;

    /*
    Construct the local FleCSolve -style Heun integrator.
    The integrator owns its temporary stage and RHS work vectors.
    */
    heun::integrator ti(
      heun::parameters(settings, flecsolve::op::ref(F), heun::make_work(U)));

    /*
    Perform:

    k1 = F(U^n)
    U* = U^n + dt k1
    k2 = F(U*)
    Uout = U^n + dt/2 (k1 + k2)
    */
    ti.advance(dt, U, Uout);

    /*
    Copy the accepted solution into HARD's conserved fields.
    */
    U.copy(Uout);

    /*
    Make HARD's primitive fields and boundary cells consistent with the final
    conserved state.
    */
    actions::flecsolve_hydro::recover_primitives_and_apply_boundaries(s, sc);
  }
};

static const auto flecsolve_advance_action = spec::register_action<control,
  state,
  hydro_flecsolve_advance,
  cp::time_integrator>();

} // namespace hard