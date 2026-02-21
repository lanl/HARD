#include "rk_stage_1.hh"

namespace hard {

template<std::size_t D>
void
RK_advance_1(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  sc.execute<tasks::rad::explicitSourceUpdate<D>>(flecsi::exec::on,
    *s.m,
    s.prim.velocity(*s.m),
    s.rad.src_t.radiation_force(*s.m),
    s.rad.src_t.radiation_pressure_tensor(*s.m),
    s.velocity_gradient(*s.m),
    //
    s.rk_dt1.total_energy_density()(*s.m),
    s.rk_dt1.momentum_energy_density()(*s.m),
    s.rad.dt_radiation_energy_density_1(*s.m));

  // We need update_u here before we compute fluxes in the presence of
  // hydro::explictSourceUpdate with body forces. See Moens'21 Eq. 24-26
  sc.execute<tasks::hydro::update_u_scalar<D>>(flecsi::exec::on,
    s.dt(*s.gt),
    //
    std::vector{
      s.rad.cons.radiation_energy_density(*s.m),
    },
    //
    std::vector{s.rad.dt_radiation_energy_density_1(*s.m)});

  using limiter = spec::limiters::weno5z;

  for(std::size_t axis = 0; axis < D; axis++) {
    sc.execute<tasks::hydro::reconstruct_primitives_scalar<D, limiter>>(
      flecsi::exec::on,
      axis,
      *s.m,
      std::vector{std::make_tuple(
        s.rad.cons.radiation_energy_density(*s.m), s.rad.f.EradFace(s.m))});

    // Calculate K1 and save it to dt_U
    sc.execute<tasks::rad::compute_interface_fluxes<D>>(flecsi::exec::on,
      axis,
      *s.m,
      s.f.hydro.uFace(s.m),
      s.f.hydro.cFace(s.m),
      s.rad.f.EradFace(s.m),
      // Riemann Fluxes
      s.rad.rf.EradF(*s.m),
      //
      s.rad.dt_radiation_energy_density_1(*s.m));
  }
}

template<std::size_t D>
void
update_vars(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  // K2 calculation in the next RK advance
  sc.execute<tasks::hydro::update_u_stage<D>>(flecsi::exec::on,
    s.dt(*s.gt),
    std::vector{s.rk_n.mass_density()(*s.m),
      s.rk_n.total_energy_density()(*s.m),
      s.rad.dt_radiation_energy_density_n(*s.m)},
    std::vector{s.rk_n.momentum_energy_density()(*s.m)},
    std::vector{s.rk_dt1.mass_density()(*s.m),
      s.rk_dt1.total_energy_density()(*s.m),
      s.rad.dt_radiation_energy_density_1(*s.m)},
    std::vector{s.rk_dt1.momentum_energy_density()(*s.m)},
    std::vector{s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.rad.cons.radiation_energy_density(*s.m)},
    std::vector{s.cons.hydro.momentum_density(*s.m)});

  // Perform primitive recovery
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

} // update_vars

template void RK_advance_1(control_policy<state, 1> &);
template void RK_advance_1(control_policy<state, 2> &);
template void RK_advance_1(control_policy<state, 3> &);

template void update_vars(control_policy<state, 1> &);
template void update_vars(control_policy<state, 2> &);
template void update_vars(control_policy<state, 3> &);

} // namespace hard
