#ifndef HARD_APPS_ACTIONS_RK_STAGE_1_HH
#define HARD_APPS_ACTIONS_RK_STAGE_1_HH

#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/interface_fluxes.hh"
#include "hydro/tasks/reconstruct.hh"
#include "hydro/tasks/rhs.hh"
#include "reconstruct_dispatch.hh"
#include "spec/reconstruct.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard::actions {

struct rk_stage_1_advance_region : annotation_action_region {
  inline static const std::string name{"rk_stage_1_advance"};
};
struct rk_stage_1_update_region : annotation_action_region {
  inline static const std::string name{"rk_stage_1_update"};
};

template<std::size_t D>
using field_double =
  flecsi::field<double>::Reference<spec::mesh<D>, spec::is::cells>;
template<std::size_t D>
using face_pair = typename hard::faces<D>::ref_pair;

template<std::size_t D, class F = std::nullptr_t>
void
advance_rk_stage_1(state<D> & s,
  flecsi::scheduler & sc,
  // conserved, dt1, face
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> v_t =
    {},
  F f = nullptr) {

  annotation_rguard<rk_stage_1_advance_region> guard;

  // RK Stage: 1 - Explicit source term (gravity) update for RT case in hydro
  // file
  annotation_rguard<tasks::external_source_task>(),
    sc.execute<tasks::external_source<D>>(flecsi::exec::on,
      *s.m,
      s.prim.velocity(*s.m),
      s.cons.hydro.mass_density(*s.m),
      s.src_t.hydro.gravity_accel(*s.m),
      // time-derivatives
      s.rk_dt1.momentum_energy_density()(*s.m),
      s.rk_dt1.total_energy_density()(*s.m));

  {
    auto scalar_v = std::vector{std::make_tuple(s.rk_dt1.mass_density()(*s.m),
                                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(s.rk_dt1.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(std::make_tuple(std::get<1>(v), std::get<0>(v)));
    }

    // We need update_u here before we compute fluxes in the presence of
    // hydro::explict_source_update with body forces. See Moens'21 Eq. 24-26
    // clang-format off
    annotation_rguard<tasks::hydro::update_u_task>(),
      sc.execute<tasks::hydro::update_u<D>>(flecsi::exec::on,
        s.dt(*s.gt),
        scalar_v,
        std::vector{
          std::make_tuple(
            s.rk_dt1.momentum_energy_density()(*s.m),
            s.cons.hydro.momentum_density(*s.m))});
    // clang-format on
  }

  // Perform primitive recovery
  annotation_rguard<tasks::hydro::conserved_to_primitive_task>(),
    sc.execute<tasks::hydro::conserved_to_primitive<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.prim.velocity(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.prim.sound_speed(*s.m),
      s.eos);

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

      dispatch_reconstruct<D>(s.reconstructor,
        sc,
        axis,
        *s.m,
        scalar_v,
        std::vector{
          std::make_tuple(s.prim.velocity(*s.m), s.f.hydro.u_face(s.m))});
    }

    annotation_rguard<tasks::hydro::get_conserved_task>(),
      sc.execute<tasks::hydro::get_conserved<D>>(flecsi::exec::on,
        *s.m,
        s.f.hydro.r_face(s.m),
        s.f.hydro.u_face(s.m),
        s.f.hydro.e_face(s.m),
        s.f.hydro.ru_face(s.m),
        s.f.hydro.re_face(s.m));

    // Calculate K1 and save it to dt_U
    annotation_rguard<tasks::hydro::compute_interface_fluxes_task>(),
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
        s.icst.gravity_const(*s.gt));

    // Apply function if present
    if constexpr(!std::is_same_v<F, std::nullptr_t>)
      f(axis);
  }
}

template<std::size_t D>
void
update_rk_stage_1(state<D> & s,
  flecsi::scheduler & sc,
  // conserved, dt1, n
  std::vector<std::tuple<field_double<D>, field_double<D>, field_double<D>>>
    v_t = {}) {

  annotation_rguard<rk_stage_1_update_region> guard;

  {
    auto scalar_v = std::vector{std::make_tuple(s.rk_n.mass_density()(*s.m),
                                  s.rk_dt1.mass_density()(*s.m),
                                  s.cons.hydro.mass_density(*s.m)),
      std::make_tuple(s.rk_n.total_energy_density()(*s.m),
        s.rk_dt1.total_energy_density()(*s.m),
        s.cons.hydro.total_energy_density(*s.m))};

    for(auto & v : v_t) {
      scalar_v.push_back(
        std::make_tuple(std::get<2>(v), std::get<1>(v), std::get<0>(v)));
    }

    // K2 calculation in the next RK advance
    // clang-format off
    annotation_rguard<tasks::hydro::update_u_stage_task>(),
      sc.execute<tasks::hydro::update_u_stage<D>>(flecsi::exec::on,
        s.dt(*s.gt),
        scalar_v,
        std::vector{
          std::make_tuple(
            s.rk_n.momentum_energy_density()(*s.m),
            s.rk_dt1.momentum_energy_density()(*s.m),
            s.cons.hydro.momentum_density(*s.m))});
    // clang-format on
  }

  // Perform primitive recovery
  annotation_rguard<tasks::hydro::conserved_to_primitive_task>(),
    sc.execute<tasks::hydro::conserved_to_primitive<D>>(flecsi::exec::on,
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
      s.prim.sound_speed(*s.m), // added but not necessary for mmhydro
      s.cons.hydro.total_energy_density(*s.m)};

    for(auto & v : v_t) {
      scalar_v.push_back(std::get<0>(v));
    }

    // Update boundary cells
    apply_tagged_boundaries<D>(sc, s, tag::boundary);
  }
}

} // namespace hard::actions

#endif //  HARD_APPS_ACTIONS_RK_STAGE_1_HH
