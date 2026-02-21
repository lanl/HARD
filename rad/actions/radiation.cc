#include "radiation.hh"

namespace hard {

template<std::size_t D>
void
radiation(control_policy<state, D> & cp) {

  using namespace flecsi;
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  sc.execute<task::rad_root::update_energy_density<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.temperature(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.rad.cons.radiation_energy_density(*s.m),
    s.rad.icst.kappa(*s.gt),
    s.dt_weighted(*s.gt),
    s.eos);

  sc.execute<tasks::rad::getGradE<D>>(flecsi::exec::on,
    *s.m,
    s.rad.cons.radiation_energy_density(*s.m),
    s.rad.limiter.gradient_rad_energy(*s.m));

  // Adaptive FLD Radiation Advance

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

  sc.execute<tasks::apply_boundaries_scalar<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.rad.limiter.lambda_bridge(*s.m)});

  sc.execute<tasks::rad::getDiff<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.rad.limiter.lambda_bridge(*s.m),
    s.rad.mgr.Diff(*s.m),
    s.rad.icst.kappa(*s.gt));

  // Initialize the diffusion coefficient
  sc.execute<tasks::rad::diffusion_init<D>>(flecsi::exec::on,
    *s.m,
    s.rad.mgr.Diff(*s.m),
    s.rad.mgr.Df_x(*s.m),
    s.rad.mgr.Df_y(*s.m),
    s.rad.mgr.Df_z(*s.m));

  // Initialize the stencil
  sc.execute<tasks::rad::stencil_init<D>>(flecsi::exec::on,
    *s.m,
    s.rad.mgr.Df_x(*s.m),
    s.rad.mgr.Df_y(*s.m),
    s.rad.mgr.Df_z(*s.m),
    s.rad.mgr.Ew(*s.m),
    s.dt_weighted(*s.gt));

  // Initialize fields
  sc.execute<tasks::rad::copy_field<D>>(flecsi::exec::on,
    *s.m,
    s.rad.cons.radiation_energy_density(*s.m),
    s.rad.mgr.Ef(*s.m));
  sc.execute<tasks::rad::const_init<D>>(
    flecsi::exec::on, *s.m, s.rad.mgr.Esf(*s.m, 1), 0.0);
  sc.execute<tasks::rad::const_init<D>>(
    flecsi::exec::on, *s.m, s.rad.mgr.Resf(*s.m), 0.0);

  std::chrono::time_point<std::chrono::system_clock> start_timer_rad =
    std::chrono::system_clock::now();

  hard::linsolve<D>(cp);

  std::chrono::time_point<std::chrono::system_clock> stop_timer_rad =
    std::chrono::system_clock::now();

  flog(info) << " Radiation Timing: "
             << (stop_timer_rad - start_timer_rad).count() * 1e-9 << " [s] "
             << std::endl;

  // Move solution from rad solver
  sc.execute<tasks::rad::copy_field<D>>(flecsi::exec::on,
    *s.m,
    s.rad.mgr.Uf(*s.m),
    s.rad.cons.radiation_energy_density(*s.m));

  // Perform primitive recovery, since energy densities have changed
  sc.execute<tasks::hydro::conservative_to_primitive<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.cons.hydro.momentum_density(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.pressure(*s.m),
    s.prim.specific_internal_energy(*s.m),
    s.prim.sound_speed(*s.m),
    s.eos);

  // Update boundary cells
  sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.rad.cons.radiation_energy_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m)},
    std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});

} // radiation_advance

// Explicit instantiation
template void radiation(control_policy<state, 1> &);
template void radiation(control_policy<state, 2> &);
template void radiation(control_policy<state, 3> &);

} // namespace hard
