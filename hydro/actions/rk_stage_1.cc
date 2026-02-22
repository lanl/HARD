#include "rk_stage_1.hh"

namespace hard {

template<std::size_t D>
void
hydro_RK_advance_1(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  // RK Stage: 1 - Explicit source term (gravity) update for RT case in hydro
  // file
  sc.execute<tasks::external_source<D>>(flecsi::exec::on,
    *s.m,
    s.prim.velocity(*s.m),
    s.src_t.hydro.gravity_force(*s.m),
    // time-derivatives
    s.rk_dt1.momentum_energy_density()(*s.m),
    s.rk_dt1.total_energy_density()(*s.m));

  // We need update_u here before we compute fluxes in the presence of
  // hydro::explictSourceUpdate with body forces. See Moens'21 Eq. 24-26
  sc.execute<tasks::hydro::update_u<D>>(flecsi::exec::on,
    s.dt(*s.gt),
    //
    std::vector{std::make_tuple(s.rk_dt1.mass_density()(*s.m),
                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(s.rk_dt1.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))},
    std::vector{std::make_tuple(s.rk_dt1.momentum_energy_density()(*s.m),
      s.cons.hydro.momentum_density(*s.m))});

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

  using limiter = spec::limiters::weno5z;

  for(std::size_t axis = 0; axis < D; axis++) {
    sc.execute<tasks::hydro::reconstruct_primitives<D, limiter>>(
      flecsi::exec::on,
      axis,
      *s.m,
      std::vector{
        std::make_tuple(s.cons.hydro.mass_density(*s.m), s.f.hydro.rFace(s.m)),
        std::make_tuple(
          s.prim.specific_internal_energy(*s.m), s.f.hydro.eFace(s.m)),
        std::make_tuple(s.prim.sound_speed(*s.m), s.f.hydro.cFace(s.m)),
        std::make_tuple(s.prim.pressure(*s.m), s.f.hydro.pFace(s.m))},
      std::vector{
        std::make_tuple(s.prim.velocity(*s.m), s.f.hydro.uFace(s.m))});

    sc.execute<tasks::hydro::reconstruct_conservatives<D>>(flecsi::exec::on,
      *s.m,
      s.f.hydro.rFace(s.m),
      s.f.hydro.uFace(s.m),
      s.f.hydro.eFace(s.m),
      s.f.hydro.ruFace(s.m),
      s.f.hydro.rEFace(s.m));

    // Calculate K1 and save it to dt_U
    sc.execute<tasks::hydro::compute_interface_fluxes<D>>(flecsi::exec::on,
      axis,
      *s.m,
      s.f.hydro.rFace(s.m),
      s.f.hydro.uFace(s.m),
      s.f.hydro.pFace(s.m),
      s.f.hydro.cFace(s.m),
      s.f.hydro.ruFace(s.m),
      s.f.hydro.rEFace(s.m),
      // Riemann Fluxes
      s.rf.hydro.rF(*s.m),
      s.rf.hydro.ruF(*s.m),
      s.rf.hydro.rEF(*s.m),

      s.rk_dt1(s.m),
      s.icst.gravity_acc(*s.gt));
  }
}

template<std::size_t D>
void
hydro_update_vars(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  // K2 calculation in the next RK advance
  sc.execute<tasks::hydro::update_u_stage<D>>(flecsi::exec::on,
    s.dt(*s.gt),
    std::vector{std::make_tuple(s.rk_n.mass_density()(*s.m),
                  s.rk_dt1.mass_density()(*s.m),
                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(s.rk_n.total_energy_density()(*s.m),
        s.rk_dt1.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))},
    std::vector{std::make_tuple(s.rk_n.momentum_energy_density()(*s.m),
      s.rk_dt1.momentum_energy_density()(*s.m),
      s.cons.hydro.momentum_density(*s.m))});

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
      s.cons.hydro.total_energy_density(*s.m)},
    std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});

} // update_vars

template void hydro_RK_advance_1(control_policy<state, 1> &);
template void hydro_RK_advance_1(control_policy<state, 2> &);
template void hydro_RK_advance_1(control_policy<state, 3> &);

template void hydro_update_vars(control_policy<state, 1> &);
template void hydro_update_vars(control_policy<state, 2> &);
template void hydro_update_vars(control_policy<state, 3> &);

} // namespace hard
