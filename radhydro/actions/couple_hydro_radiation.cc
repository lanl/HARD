#include "../modules/rad/tasks/utils.hh"
#include "state.hh"

#include "../modules/rad/tasks/rad.hh"

namespace hard {

template<std::size_t D>
void
couple_hydro_radiation(control_policy<state, D> & cp) {

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

} // hydro_couple_radiation

inline control<state, 1>::action<couple_hydro_radiation<1>,
  cp::couple_hydro_radiation_1>
  couple_hydro_radiation_1_1d;
inline control<state, 2>::action<couple_hydro_radiation<2>,
  cp::couple_hydro_radiation_1>
  couple_hydro_radiation_1_2d;
inline control<state, 3>::action<couple_hydro_radiation<3>,
  cp::couple_hydro_radiation_1>
  couple_hydro_radiation_1_3d;

inline control<state, 1>::action<couple_hydro_radiation<1>,
  cp::couple_hydro_radiation_2>
  couple_hydro_radiation_2_1d;
inline control<state, 2>::action<couple_hydro_radiation<2>,
  cp::couple_hydro_radiation_2>
  couple_hydro_radiation_2_2d;
inline control<state, 3>::action<couple_hydro_radiation<3>,
  cp::couple_hydro_radiation_2>
  couple_hydro_radiation_2_3d;

} // namespace hard
