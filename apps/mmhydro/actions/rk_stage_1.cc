#include "state.hh"

#include "gfm/tasks/cons2prim.hh"
#include "gfm/tasks/extrapolation.hh"
#include "gfm/tasks/interface_fluxes_gfm.hh"
#include "gfm/tasks/reconstruct.hh"
#include "gfm/tasks/utils.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/limiter.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard {

template<std::size_t D>
using field_double =
  flecsi::field<double>::Reference<spec::mesh<D>, spec::is::cells>;
template<std::size_t D>
using face_pair = typename hard::faces<D>::ref_pair;

template<std::size_t D, class F = std::nullptr_t>
void
advance_rk_stage_1(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt1, face
  std::vector<std::tuple<field_double<D>, field_double<D>, face_pair<D>>> v_t =
    {},
  F f = nullptr) {

  // RK Stage: 1 - Explicit source term (gravity) update for RT case in hydro
  // file
  sc.execute<tasks::external_source<D>>(flecsi::exec::on,
    *s.m,
    s.prim.velocity(*s.m),
    s.src_t.hydro.gravity_force(*s.m),
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
  sc.execute<tasks::gfm::conservative_to_primitive_mm<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.cons.hydro.momentum_density(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.pressure(*s.m),
    s.prim.specific_internal_energy(*s.m),
    s.prim.sound_speed(*s.m),
    s.ls.material(*s.m),
    s.eos_mm);

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

      // TODO : this is the where extrapolation should go to feel ghost

      // Compute unit normal of level set gradient
      // dp = grad(psi)
      sc.execute<tasks::extrapolation::gradient_centered<D>>(
        flecsi::exec::on, *s.m, s.ls.psi(*s.m), s.extp.dp(*s.m));

      // normal = dp/|dp|
      sc.execute<tasks::extrapolation::unit_normals<D>>(
        flecsi::exec::on, *s.m, s.extp.dp(*s.m), s.extp.normal(*s.m));

      // Stage 1: rhs(q^n) -> rhs
      sc.execute<tasks::extrapolation::rhs_constant_extrapolation_minmod<D>>(
        flecsi::exec::on,
        *s.m,
        s.extp.q_scalar(*s.m),
        s.ls.psi(*s.m),
        s.extp.normal(*s.m),
        s.extp.q_rhs(*s.m));

      // q1 = q^n + dt * rhs(q^n)
      sc.execute<tasks::extrapolation::advance_q_rk2<D, 1>>(flecsi::exec::on,
        *s.m,
        s.ls.psi(*s.m),
        s.extp.q_scalar_known(*s.m),
        s.extp.q_scalar(*s.m),
        s.extp.q_rhs(*s.m),
        s.extp.q_scalar2(*s.m), // dummy (unused in STEP=1)
        s.extp.q_scalar1(*s.m), // output stage-1
        s.dt(*s.gt)); // dtau in the origianl work

      // clang-format off
    sc.execute<tasks::hydro::reconstruct_primitives<D, limiter>>(
      flecsi::exec::on,
      axis,
      *s.m,
      s.ls.material(*s.m),
      scalar_v,
      std::vector{
        std::make_tuple(
          s.prim.velocity(*s.m),
          s.f.hydro.u_face(s.m))});
      // clang-format on
    }

    // Stage 2: rhs(q1) -> rhs1
    sc.execute<tasks::extrapolation::rhs_constant_extrapolation_minmod<D>>(
      flecsi::exec::on,
      *s.m,
      s.extp.q_scalar1(*s.m),
      s.ls.psi(*s.m),
      s.extp.normal(*s.m),
      s.extp.q_rhs1(*s.m));

    // q^{n+1} = 0.5*( q^n + q1 + dt*rhs(q1) )
    sc.execute<tasks::extrapolation::advance_q_rk2<D, 2>>(flecsi::exec::on,
      *s.m,
      s.ls.psi(*s.m),
      s.extp.q_scalar_known(*s.m),
      s.extp.q_scalar(*s.m), // q^n
      s.extp.q_rhs1(*s.m), // rhs(q1)
      s.extp.q_scalar1(*s.m), // q_stage
      s.extp.q_scalar2(*s.m), // q_new
      s.dt(*s.gt)); // dtau in original work

    // Swap
    sc.execute<hard::tasks::util::copy_scal_field<dm::quantities, D, double>>(
      flecsi::exec::on, *s.m, s.extp.q_scalar2(*s.m), s.extp.q_scalar(*s.m));

    sc.execute<tasks::hydro::reconstruct_conservatives<D>>(flecsi::exec::on,
      *s.m,
      s.f.hydro.r_face(s.m),
      s.f.hydro.u_face(s.m),
      s.f.hydro.e_face(s.m),
      s.f.hydro.ru_face(s.m),
      s.f.hydro.re_face(s.m));

    // Calculate K1 and save it to dt_U
    sc.execute<tasks::gfm::compute_interface_fluxes_gfm<D>>(flecsi::exec::on,
      axis,
      *s.m,
      s.f.hydro.r_face(s.m),
      s.f.hydro.u_face(s.m),
      s.f.hydro.p_face(s.m),
      s.f.hydro.c_face(s.m),
      s.f.hydro.ru_face(s.m),
      s.f.hydro.re_face(s.m),
      // GFM
      s.ls.material(*s.m),
      s.extp.f.normal_face(s.m),
      // Riemann Fluxes
      s.rf.hydro.r_f(*s.m),
      s.rf.hydro.ru_f(*s.m),
      s.rf.hydro.re_f(*s.m),
      s.rk_dt1(s.m),
      s.icst.gravity_acc(*s.gt),
      s.eos_mm);

    // Apply function if present
    if constexpr(!std::is_same_v<F, std::nullptr_t>)
      f(axis);
  }
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
  sc.execute<tasks::gfm::conservative_to_primitive_mm<D>>(flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.cons.hydro.momentum_density(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.pressure(*s.m),
    s.prim.specific_internal_energy(*s.m),
    s.prim.sound_speed(*s.m),
    s.ls.material(*s.m),
    s.eos_mm);

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
