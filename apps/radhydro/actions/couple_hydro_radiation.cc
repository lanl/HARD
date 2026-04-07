#include "rad/tasks/utils.hh"
#include "state.hh"

#include "rad/tasks/rad.hh"

#include <spec/runtime.hh>

namespace hard {

template<std::size_t D>
struct couple_hydro_radiation {
  static void action(control_policy<state, D> & cp) {
    auto & s = cp.state();
    flecsi::scheduler & sc = cp.scheduler();

    sc.execute<tasks::rad::getGradE<D>>(flecsi::exec::on,
      *s.m,
      s.rad.cons.radiation_energy_density(*s.m),
      s.rad.limiter.gradient_rad_energy(*s.m));

    // Standard (Constant) FLD formulation

    sc.execute<tasks::rad::getLambda<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.rad.limiter.gradient_rad_energy(*s.m),
      s.rad.limiter.magnitude_gradient_rad_energy(*s.m),
      s.rad.limiter.R_value(*s.m),
      s.rad.limiter.lambda_bridge(*s.m),
      s.rad.icst.kappa(*s.gt),
      s.rad.icst.limiter_id(*s.gt));

    sc.execute<tasks::rad::getEddFactor<D>>(flecsi::exec::on,
      *s.m,
      s.rad.limiter.lambda_bridge(*s.m),
      s.rad.limiter.eddington_factor(*s.m),
      s.rad.icst.limiter_id(*s.gt),
      s.rad.icst.closure_id(*s.gt));

    sc.execute<tasks::rad::getTensorP<D>>(flecsi::exec::on,
      *s.m,
      s.rad.src_t.radiation_pressure_tensor(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.rad.limiter.gradient_rad_energy(*s.m),
      s.rad.limiter.magnitude_gradient_rad_energy(*s.m),
      s.rad.limiter.eddington_factor(*s.m));

    sc.execute<tasks::rad::getRadForce<D>>(flecsi::exec::on,
      *s.m,
      s.rad.limiter.lambda_bridge(*s.m),
      s.rad.limiter.gradient_rad_energy(*s.m),
      s.rad.src_t.radiation_force(*s.m));

    sc.execute<tasks::rad::getGradV<D>>(
      flecsi::exec::on, *s.m, s.velocity_gradient(*s.m), s.prim.velocity(*s.m));
  }
};

static const auto couple_1_action = spec::register_action<control,
  state,
  couple_hydro_radiation,
  cp::couple_hydro_radiation_1>();

static const auto couple_2_action = spec::register_action<control,
  state,
  couple_hydro_radiation,
  cp::couple_hydro_radiation_2>();

} // namespace hard
