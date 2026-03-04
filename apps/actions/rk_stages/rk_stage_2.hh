#ifndef HARD_APPS_ACTIONS_RK_STAGE_2_HH
#define HARD_APPS_ACTIONS_RK_STAGE_2_HH

#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/interface_fluxes.hh"
#include "hydro/tasks/reconstruct.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard::actions {

template<std::size_t D>
using field_double =
  flecsi::field<double>::Reference<spec::mesh<D>, spec::is::cells>;
template<std::size_t D>
using face_pair = typename hard::faces<D>::ref_pair;

template<std::size_t D, class F = std::nullptr_t>
void
advance_rk_stage_2(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt2, face
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> v_t =
    {},
  F f = nullptr) {

  // RK Stage: 2 - Explicit source term (gravity) update for RT case in hydro
  // file
  sc.execute<tasks::external_source<D>>(flecsi::exec::on,
    *s.m,
    s.prim.velocity(*s.m),
    s.src_t.hydro.gravity_force(*s.m),
    // time-derivatives
    s.rk_dt2.momentum_energy_density()(*s.m),
    s.rk_dt2.total_energy_density()(*s.m));

  {
    auto scalar_v = std::vector{std::make_tuple(s.rk_dt2.mass_density()(*s.m),
                                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(s.rk_dt2.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(std::make_tuple(std::get<1>(v), std::get<0>(v)));
    }

    // We need update_u here before we compute fluxes in the presence of
    // hydro::explict_source_update with body forces. See Moens'21 Eq. 24-26
    sc.execute<tasks::hydro::update_u<D>>(flecsi::exec::on,
      s.dt(*s.gt),
      scalar_v,
      std::vector{std::make_tuple(s.rk_dt2.momentum_energy_density()(*s.m),
        s.cons.hydro.momentum_density(*s.m))});
  }

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
    {
      auto scalar_v = std::vector{
        std::make_tuple(s.cons.hydro.mass_density(*s.m), s.f.hydro.r_face(s.m)),
        std::make_tuple(
          s.prim.specific_internal_energy(*s.m), s.f.hydro.e_face(s.m)),
        std::make_tuple(s.prim.sound_speed(*s.m), s.f.hydro.c_face(s.m)),
        std::make_tuple(s.prim.pressure(*s.m), s.f.hydro.p_face(s.m))};

      for(auto & v : v_t) {
        scalar_v.push_back(std::make_tuple(std::get<0>(v), std::get<2>(v)));
      }

      // clang-format off
    sc.execute<tasks::hydro::reconstruct_primitives<D, limiter>>(
      flecsi::exec::on,
      axis,
      *s.m,
      scalar_v,
      std::vector{
        std::make_tuple(
          s.prim.velocity(*s.m),
          s.f.hydro.u_face(s.m))});
      // clang-format on
    }

    sc.execute<tasks::hydro::reconstruct_conservatives<D>>(flecsi::exec::on,
      *s.m,
      s.f.hydro.r_face(s.m),
      s.f.hydro.u_face(s.m),
      s.f.hydro.e_face(s.m),
      s.f.hydro.ru_face(s.m),
      s.f.hydro.re_face(s.m));

    // Calculate K1 and save it to dt_U
    sc.execute<tasks::hydro::compute_interface_fluxes<D>>(flecsi::exec::on,
      axis,
      *s.m,
      s.f.hydro.r_face(s.m),
      s.f.hydro.u_face(s.m),
      s.f.hydro.p_face(s.m),
      s.f.hydro.c_face(s.m),
      s.f.hydro.ru_face(s.m),
      s.f.hydro.re_face(s.m),
      // Riemann Fluxes
      s.rf.hydro.r_f(*s.m),
      s.rf.hydro.ru_f(*s.m),
      s.rf.hydro.re_f(*s.m),
      s.rk_dt1(s.m),
      s.icst.gravity_acc(*s.gt));

    // Apply function if present
    if constexpr(!std::is_same_v<F, std::nullptr_t>)
      f(axis);
  }
}

template<std::size_t D>

void
update_rk_stage_2(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt1, dt2, n
  std::vector<std::
      tuple<field_double<D>, field_double<D>, field_double<D>, field_double<D>>>
    v_t = {}) {

  {
    auto scalar_v = std::vector{std::make_tuple(s.rk_dt1.mass_density()(*s.m),
                                  s.rk_dt2.mass_density()(*s.m)),
      std::make_tuple(s.rk_dt1.total_energy_density()(*s.m),
        s.rk_dt2.total_energy_density()(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(std::make_tuple(std::get<1>(v), std::get<2>(v)));
    }

    // First compute K1' = (K1 + K2) * 0.5
    // clang-format off
  sc.execute<tasks::hydro::add_k1_k2<D>>(flecsi::exec::on,
    scalar_v,
    std::vector{
      std::make_tuple(
        s.rk_dt1.momentum_energy_density()(*s.m),
        s.rk_dt2.momentum_energy_density()(*s.m))});
    // clang-format on
  }

  {
    auto scalar_v = std::vector{std::make_tuple(s.rk_dt1.mass_density()(*s.m),
                                  s.rk_n.mass_density()(*s.m)),
      std::make_tuple(s.rk_dt1.total_energy_density()(*s.m),
        s.rk_n.total_energy_density()(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(std::make_tuple(std::get<1>(v), std::get<3>(v)));
    }

    // Now get U_n(+1) = U_n + h * K1'
    sc.execute<tasks::hydro::update_u<D>>(flecsi::exec::on,
      s.dt(*s.gt),
      scalar_v,
      std::vector{std::make_tuple(s.rk_dt1.momentum_energy_density()(*s.m),
        s.rk_n.momentum_energy_density()(*s.m))});
  }

  {
    auto scalar_v = std::vector{std::make_tuple(s.rk_n.mass_density()(*s.m),
                                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(s.rk_n.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(std::make_tuple(std::get<3>(v), std::get<0>(v)));
    }

    // Finish by updating the values stored in U_n to U
    sc.execute<tasks::hydro::store_current_state<D>>(flecsi::exec::on,
      scalar_v,
      std::vector{std::make_tuple(s.rk_n.momentum_energy_density()(*s.m),
        s.cons.hydro.momentum_density(*s.m))});
  }

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

  {
    auto scalar_v = std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.cons.hydro.total_energy_density(*s.m)};

    for(auto & v : v_t) {
      scalar_v.push_back(std::get<0>(v));
    }

    // Update boundary cells
    sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
      *s.m,
      s.icst.bmap(*s.gt),
      scalar_v,
      std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});
  }
}

} // namespace hard::actions

#endif //  HARD_APPS_ACTIONS_RK_STAGE_2_HH
