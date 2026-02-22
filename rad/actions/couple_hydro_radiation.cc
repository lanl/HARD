#include "couple_hydro_radiation.hh"

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

template void couple_hydro_radiation(control_policy<state, 1> &);
template void couple_hydro_radiation(control_policy<state, 2> &);
template void couple_hydro_radiation(control_policy<state, 3> &);

} // namespace hard
