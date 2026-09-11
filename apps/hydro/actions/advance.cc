#include "state.hh"

#include "../../actions/flecsolve/hydro_flecsolve.hh"
#include "../../actions/flecsolve/solver.hh"
#include "flecsolve/time-integrators/rk23.hh"
#include "flecsolve/time-integrators/rk45.hh"

#include <flecsi/flog.hh>
#include <spec/runtime.hh>

namespace hard {

struct hydro_advance_region : annotation_action_region {
  inline static const std::string name{"hydro_advance"};
};

namespace solver = time_integrator::solver;
namespace rk23 = flecsolve::time_integrator::rk23;
namespace rk45 = flecsolve::time_integrator::rk45;

template<std::size_t D>
struct hydro_advance {
  static void action(control_policy<state, D> & cp) {

    annotation_rguard<hydro_advance_region> guard;

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
    solver step.
    */
    auto U = actions::flecsolve_hydro::make_state_vector(s);
    auto Uout = actions::flecsolve_hydro::make_saved_state_vector(s);

    /*
    Configure one fixed solver step over the timestep selected by HARD
    */

    const bool flecsolve_solver =
      (s.use_solver == spec::slv::rk23 || s.use_solver == spec::slv::rk45);

    // Solver (time integrator)
    if(flecsolve_solver) {
      if(s.use_solver == spec::slv::rk23) {
        rk23::settings settings{};

        settings.initial_time = 0.0;
        settings.final_time = dt;
        settings.max_steps = 1;

        settings.initial_dt = dt;
        settings.min_dt = dt;
        settings.max_dt = dt;
        settings.use_fixed_dt = true;

        settings.safety_factor = 0.9f;
        settings.atol = 1.0e-12f;

        rk23::integrator ti(rk23::parameters(
          settings, flecsolve::op::ref(F), rk23::make_work(U)));

        ti.advance(dt, U, Uout);
      }
      else if(s.use_solver == spec::slv::rk45) {
        rk45::settings settings{};

        settings.initial_time = 0.0;
        settings.final_time = dt;
        settings.max_steps = 1;

        settings.initial_dt = dt;
        settings.min_dt = dt;
        settings.max_dt = dt;
        settings.use_fixed_dt = true;

        settings.safety_factor = 0.9f;
        settings.atol = 1.0e-12f;

        rk45::integrator ti(rk45::parameters(
          settings, flecsolve::op::ref(F), rk45::make_work(U)));

        ti.advance(dt, U, Uout);
      }
    }
    else {

      solver::settings settings{};
      settings.initial_time = 0.0;
      settings.final_time = dt;
      settings.max_steps = 1;

      settings.initial_dt = dt;
      settings.min_dt = dt;
      settings.max_dt = dt;

      if(s.use_solver == spec::slv::euler) {
        auto ti =
          solver::make_integrator<spec::slv::euler>(solver::parameters(settings,
            flecsolve::op::ref(F),
            solver::make_work<spec::slv::euler_vecs>(U)));

        ti.advance(dt, U, Uout);
      }
      else if(s.use_solver == spec::slv::heun) {
        auto ti =
          solver::make_integrator<spec::slv::heun>(solver::parameters(settings,
            flecsolve::op::ref(F),
            solver::make_work<spec::slv::heun_vecs>(U)));

        ti.advance(dt, U, Uout);
      }
      else if(s.use_solver == spec::slv::rk4) {
        auto ti =
          solver::make_integrator<spec::slv::rk4>(solver::parameters(settings,
            flecsolve::op::ref(F),
            solver::make_work<spec::slv::rk4_vecs>(U)));

        ti.advance(dt, U, Uout);
      }
    }

    // Copy the accepted solution into HARD's conserved fields.
    U.copy(Uout);

    // Update primitives and boundaries
    actions::flecsolve_hydro::recover_primitives_and_apply_boundaries(s, sc);
  }
};

static const auto advance_action =
  spec::register_action<control, state, hydro_advance, cp::time_integrator>();

} // namespace hard
