#include "rk_stages.hh"

namespace hard {

template<std::size_t D, class F = std::nullptr_t>
void
advance_rk_stage_1(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt1, face
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> v_t =
    {},
  F f = nullptr) {

  // Stage 1 uses rk_dt1 for the explicit source/update and recomputes the
  // sound speed from EOS in ghost cells.
  advance_rk_stage_common(s, sc, s.rk_dt1, /*do_sound_speed=*/true, v_t, f);
}

template<std::size_t D>
void
update_rk_stage_1(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt1, n
  std::vector<std::tuple<field_double<D>, field_double<D>, field_double<D>>>
    v_t = {}) {

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
hydro_RK_advance_1(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  advance_rk_stage_1(s, sc);
}

template<std::size_t D>
void
hydro_update_vars(control_policy<state, D> & cp) {

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

  update_rk_stage_1(s, sc);

} // update_vars

inline control<state, 1>::action<hydro_RK_advance_1<1>, cp::rk_stage_1>
  hydro_rk_stage_1_1d;
inline control<state, 2>::action<hydro_RK_advance_1<2>, cp::rk_stage_1>
  hydro_rk_stage_1_2d;
inline control<state, 3>::action<hydro_RK_advance_1<3>, cp::rk_stage_1>
  hydro_rk_stage_1_3d;

inline control<state, 1>::action<hydro_update_vars<1>, cp::rk_stage_1>
  hydro_rk_stage_1_update_1d;
inline control<state, 2>::action<hydro_update_vars<2>, cp::rk_stage_1>
  hydro_rk_stage_1_update_2d;
inline control<state, 3>::action<hydro_update_vars<3>, cp::rk_stage_1>
  hydro_rk_stage_1_update_3d;

inline const auto hydro_dep_update_1d =
  hydro_rk_stage_1_update_1d.add(hydro_rk_stage_1_1d);
inline const auto hydro_dep_update_2d =
  hydro_rk_stage_1_update_2d.add(hydro_rk_stage_1_2d);
inline const auto hydro_dep_update_3d =
  hydro_rk_stage_1_update_3d.add(hydro_rk_stage_1_3d);

} // namespace hard
