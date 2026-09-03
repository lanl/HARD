#include "state.hh"

#include "../../actions/flecsolve/hydro_flecsolve.hh"

#include <flecsi/flog.hh>
#include <spec/runtime.hh>

#include "flecsolve/time-integrators/rk23.hh"

namespace hard {

struct hydro_flecsolve_advance_region : annotation_action_region {
  inline static const std::string name{"hydro_flecsolve_advance_region"};
};

namespace rk23 = flecsolve::time_integrator::rk23;

template<std::size_t D>
struct hydro_flecsolve_advance {
  static void action(control_policy<state, D> & cp) {

    annotation_rguard<hydro_flecsolve_advance_region> guard;

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
      NOTE:

      The rhs_operator now explicitly specifies its input/output variables
      (mass_density, momentum_density, total_energy_density), so wrapping
      with flecsolve::op::core should work correctly. However, we use the
      operator directly here without wrapping.
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

    rk23::integrator ti(
      rk23::parameters(settings, flecsolve::op::ref(F), rk23::make_work(U)));

    ti.advance(dt, U, Uout);

    U.copy(Uout);

    actions::flecsolve_hydro::recover_primitives_and_apply_boundaries(s, sc);
  }
};

static const auto flecsolve_advance_action = spec::register_action<control,
  state,
  hydro_flecsolve_advance,
  cp::time_integrator>();

} // namespace hard
