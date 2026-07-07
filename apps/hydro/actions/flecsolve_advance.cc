#include "state.hh"

#include "../../actions/flecsolve/hydro_flecsolve.hh"

#include <spec/runtime.hh>
#include <flecsi/flog.hh>

#include "flecsolve/time-integrators/rk23.hh"

namespace hard {

namespace rk23 = flecsolve::time_integrator::rk23;

template<std::size_t D>
struct hydro_flecsolve_advance {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    const double dt =
      sc.template execute<actions::flecsolve_hydro::read_dt_task>(
          flecsi::exec::on, s.dt(*s.gt))
        .get();

    if(dt <= 0.0) {
      return;
    }

    /*
      IMPORTANT:

      Do not wrap this in flecsolve::op::core. The core wrapper tries to
      subset anonymous FleCSolve variables, but our hydro multivector uses
      hard::actions::flecsolve_hydro::var.
    */
    actions::flecsolve_hydro::rhs_params<D> params{&s, &sc};
    actions::flecsolve_hydro::rhs_operator<D> F(params);

    auto U = actions::flecsolve_hydro::make_state_vector(s);
    auto Uout = actions::flecsolve_hydro::make_saved_state_vector(s);

    rk23::settings settings{};

    settings.initial_time = 0.0;
    settings.final_time = dt;
    settings.max_steps = 1;

    settings.initial_dt = dt;
    settings.min_dt = dt;
    settings.max_dt = dt;

    settings.safety_factor = 0.9f;
    settings.atol = 1.0e-12f;
    settings.use_fixed_dt = true;

    rk23::integrator ti(rk23::parameters(
      settings,
      flecsolve::op::ref(F),
      rk23::make_work(U)));

    ti.advance(dt, U, Uout);

    U.copy(Uout);

    actions::flecsolve_hydro::recover_primitives_and_apply_boundaries(s, sc);
  }
};

static const auto flecsolve_advance_action =
  spec::register_action<control, state, hydro_flecsolve_advance, cp::rk_stage_1>();

} // namespace hard