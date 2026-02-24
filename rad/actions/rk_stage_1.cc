#include "state.hh"
#include "utils.hh"

#include "../modules/hydro/tasks/cons2prim.hh"
#include "../modules/hydro/tasks/external_source.hh"
#include "../modules/hydro/tasks/interface_fluxes.hh"
#include "../modules/hydro/tasks/reconstruct.hh"
#include "../modules/hydro/tasks/time_derivative.hh"
#include "../modules/rad/tasks/interface_fluxes.hh"
#include "../modules/rad/tasks/rad.hh"
#include "../modules/spec/limiter.hh"
#include "../modules/spec/tasks/boundaries/boundary.hh"

namespace hard {

template<std::size_t D>
void
RK_advance_1(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  sc.execute<tasks::rad::explicit_source_update<D>>(flecsi::exec::on,
    *s.m,
    s.prim.velocity(*s.m),
    s.rad.src_t.radiation_force(*s.m),
    s.rad.src_t.radiation_pressure_tensor(*s.m),
    s.velocity_gradient(*s.m),
    //
    s.rk_dt1.total_energy_density()(*s.m),
    s.rk_dt1.momentum_energy_density()(*s.m),
    s.rad.dt_radiation_energy_density_1(*s.m));

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
  // hydro::explict_source_update with body forces. See Moens'21 Eq. 24-26
  // clang-format off
  sc.execute<tasks::hydro::update_u<D>>(flecsi::exec::on,
    s.dt(*s.gt),
    std::vector{
      std::make_tuple(
        s.rk_dt1.mass_density()(*s.m),
        s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(
        s.rk_dt1.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m)),
      std::make_tuple(
        s.rad.dt_radiation_energy_density_1(*s.m),
        s.rad.cons.radiation_energy_density(*s.m))},
    std::vector{
      std::make_tuple(
        s.rk_dt1.momentum_energy_density()(*s.m),
        s.cons.hydro.momentum_density(*s.m))});
  // clang-format on

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
    // clang-format off
    sc.execute<tasks::hydro::reconstruct_primitives<D, limiter>>(
      flecsi::exec::on,
      axis,
      *s.m,
      std::vector{
        std::make_tuple(
          s.cons.hydro.mass_density(*s.m),
          s.f.hydro.rFace(s.m)),
        std::make_tuple(
          s.prim.specific_internal_energy(*s.m),
          s.f.hydro.eFace(s.m)),
        std::make_tuple(
          s.prim.sound_speed(*s.m),
          s.f.hydro.cFace(s.m)),
        std::make_tuple(
          s.rad.cons.radiation_energy_density(*s.m),
          s.rad.f.EradFace(s.m)),
        std::make_tuple(
          s.prim.pressure(*s.m),
          s.f.hydro.pFace(s.m))},
      std::vector{
        std::make_tuple(
          s.prim.velocity(*s.m),
          s.f.hydro.uFace(s.m))});
    // clang-format on

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
  // clang-format off
  sc.execute<tasks::hydro::update_u_stage_scalar<D>>(flecsi::exec::on,
    s.dt(*s.gt),
    std::vector{
      std::make_tuple(
        s.rad.dt_radiation_energy_density_n(*s.m),
        s.rad.dt_radiation_energy_density_1(*s.m),
        s.rad.cons.radiation_energy_density(*s.m))});
  // clang-format on

  // Update boundary cells
  sc.execute<tasks::apply_boundaries_scalar<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.rad.cons.radiation_energy_density(*s.m)});

} // update_vars

inline control<state, 1>::action<RK_advance_1<1>, cp::rk_stage_1>
  rad_rk_stage_1_1d;
inline control<state, 2>::action<RK_advance_1<2>, cp::rk_stage_1>
  rad_rk_stage_1_2d;
inline control<state, 3>::action<RK_advance_1<3>, cp::rk_stage_1>
  rad_rk_stage_1_3d;

inline control<state, 1>::action<update_vars<1>, cp::rk_stage_1>
  rad_rk_stage_1_update_1d;
inline control<state, 2>::action<update_vars<2>, cp::rk_stage_1>
  rad_rk_stage_1_update_2d;
inline control<state, 3>::action<update_vars<3>, cp::rk_stage_1>
  rad_rk_stage_1_update_3d;

//  hydro_rk_stage_1 -> rad_rk_stage_1
inline const auto dep_update_1d =
  rad_rk_stage_1_update_1d.add(rad_rk_stage_1_1d);
inline const auto dep_update_2d =
  rad_rk_stage_1_update_2d.add(rad_rk_stage_1_2d);
inline const auto dep_update_3d =
  rad_rk_stage_1_update_3d.add(rad_rk_stage_1_3d);

} // namespace hard
