#include "rk_stages.hh"

namespace hard {

template<std::size_t D, class F = std::nullptr_t>
void
advance_rk_stage_2(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt2, face
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> v_t =
    {},
  F f = nullptr) {

  // Stage 2 uses rk_dt2 for the explicit source/update and does NOT recompute
  // the sound speed from EOS in ghost cells.
  advance_rk_stage_common(s, sc, s.rk_dt2, /*do_sound_speed=*/false, v_t, f);
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
    sc.execute<common::tasks::utils::copy_scalar_v_vector<D, vec<D>>>(
      flecsi::exec::on,
      scalar_v,
      std::vector{std::make_tuple(s.rk_n.momentum_energy_density()(*s.m),
        s.cons.hydro.momentum_density(*s.m))});
  }

  // Perform primitive recovery
  perform_cons2prim(s, sc);

  {
    auto scalar_v = std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.cons.hydro.total_energy_density(*s.m)};

    for(auto & v : v_t) {
      scalar_v.push_back(std::get<0>(v));
    }

    // Update material/psi ghost cells and then boundary cells.
    apply_scalar_boundaries(s, sc, scalar_v);
  }
}

template<std::size_t D>
void
RK_advance_2(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  advance_rk_stage_2(s, sc);
}

template<std::size_t D>
void
update_vars_2(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  update_rk_stage_2(s, sc);

} // update_vars

inline control<state, 1>::action<RK_advance_2<1>, cp::rk_stage_2>
  hydro_rk_stage_2_1d;
inline control<state, 2>::action<RK_advance_2<2>, cp::rk_stage_2>
  hydro_rk_stage_2_2d;
inline control<state, 3>::action<RK_advance_2<3>, cp::rk_stage_2>
  hydro_rk_stage_2_3d;

inline control<state, 1>::action<update_vars_2<1>, cp::rk_stage_2>
  hydro_rk_stage_2_update_1d;
inline control<state, 2>::action<update_vars_2<2>, cp::rk_stage_2>
  hydro_rk_stage_2_update_2d;
inline control<state, 3>::action<update_vars_2<3>, cp::rk_stage_2>
  hydro_rk_stage_2_update_3d;

inline const auto dep_update_2_1d =
  hydro_rk_stage_2_update_1d.add(hydro_rk_stage_2_1d);
inline const auto dep_update_2_2d =
  hydro_rk_stage_2_update_2d.add(hydro_rk_stage_2_2d);
inline const auto dep_update_2_3d =
  hydro_rk_stage_2_update_3d.add(hydro_rk_stage_2_3d);

} // namespace hard
