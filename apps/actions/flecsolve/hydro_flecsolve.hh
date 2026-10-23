#ifndef HARD_APPS_ACTIONS_FLECSOLVE_HYDRO_FLECSOLVE_HH
#define HARD_APPS_ACTIONS_FLECSOLVE_HYDRO_FLECSOLVE_HH

#include "flecsolve/operators/core.hh"
#include "flecsolve/vectors/multi.hh"
#include "flecsolve/vectors/topo_view.hh"

#include "../rk_stages/reconstruct_dispatch.hh"
#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/external_source.hh"
#include "hydro/tasks/interface_fluxes.hh"
#include "hydro/tasks/reconstruct.hh"
#include "spec/tasks/boundaries/boundary.hh"

namespace hard::actions::flecsolve_hydro {

enum class var { mass_density, total_energy_density, momentum_density };

struct make_conserved_vector_task : annotation_action_region {
  inline static const std::string name{"make conserved vector"};
};

// Helper function to make the repeated conserved vector below
FLECSI_INLINE_TARGET auto
make_conserved_vector(auto rho, auto e, auto mom) {

  annotation_rguard<make_conserved_vector_task> guard;

  return flecsolve::vec::make(
    flecsolve::vec::make(flecsolve::variable<var::mass_density>, rho),
    flecsolve::vec::make(flecsolve::variable<var::total_energy_density>, e),
    flecsolve::vec::make(flecsolve::variable<var::momentum_density>, mom));
}

struct make_state_vector_task : annotation_action_region {
  inline static const std::string name{"make_state_vector"};
};

/*
  FleCSolve view of the main HARD hydro conserved state:

      U = (rho, rho*u, rho*E)

  These are views into s.cons.hydro.
*/
template<std::size_t D>
auto
make_state_vector(state<D> & s) {
  auto & m = *s.m;

  annotation_rguard<make_state_vector_task> guard;

  return make_conserved_vector(s.cons.hydro.mass_density(m),
    s.cons.hydro.total_energy_density(m),
    s.cons.hydro.momentum_density(m));
}

struct make_rhs_vector_task : annotation_action_region {
  inline static const std::string name{"make_rhs_vector"};
};

/*
  Temporary RHS vector.

  For now, we use s.rk_dt1 as FleCSolve RHS storage:

      R = dU/dt = (d rho/dt, d rho*u/dt, d rho*E/dt)
*/
template<std::size_t D>
auto
make_rhs_vector(state<D> & s) {
  auto & m = *s.m;

  annotation_rguard<make_rhs_vector_task> guard;

  return make_conserved_vector(s.rk_dt1.mass_density()(m),
    s.rk_dt1.total_energy_density()(m),
    s.rk_dt1.momentum_energy_density()(m));
}

struct make_rhs2_vector_task : annotation_action_region {
  inline static const std::string name{"make_rhs2_vector"};
};

/*
  Temporary second-stage RHS vector.

  This uses s.rk_dt2 as the stage-2 RHS storage:

      R2 = F(U*)
*/
template<std::size_t D>
auto
make_rhs2_vector(state<D> & s) {
  auto & m = *s.m;

  annotation_rguard<make_rhs2_vector_task> guard;

  return make_conserved_vector(s.rk_dt2.mass_density()(m),
    s.rk_dt2.total_energy_density()(m),
    s.rk_dt2.momentum_energy_density()(m));
}

struct make_saved_vector_task : annotation_action_region {
  inline static const std::string name{"make_saved_vector"};
};

/*
  FleCSolve view of the saved old state U^n.

  HARD's rhs action stores the beginning-of-step state in s.rk_n before
  the RK update starts.
*/
template<std::size_t D>
auto
make_saved_state_vector(state<D> & s) {
  auto & m = *s.m;

  annotation_rguard<make_saved_vector_task> guard;

  return make_conserved_vector(s.rk_n.mass_density()(m),
    s.rk_n.total_energy_density()(m),
    s.rk_n.momentum_energy_density()(m));
}

// Return the dt future
inline double
read_dt_task(flecsi::exec::cpu,
  single<double>::accessor<flecsi::ro> dt) noexcept {
  return *dt;
}

/*
  After we update the conserved state U, we need to recompute primitive
  variables and refresh boundary cells, just like HARD's old update stages did.
*/
template<std::size_t D>
void
recover_primitives_and_apply_boundaries(state<D> & s, flecsi::scheduler & sc) {
  annotation_rguard<tasks::hydro::conserved_to_primitive_task>(),
    sc.template execute<tasks::hydro::conserved_to_primitive<D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.prim.velocity(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.prim.sound_speed(*s.m),
      s.eos);

  apply_tagged_boundaries<D>(sc, s, tag::boundary);
}

template<std::size_t D>
struct rhs_params {
  state<D> * s;
  flecsi::scheduler * sc;
};

template<std::size_t D>
struct rhs_operator : flecsolve::op::base<rhs_params<D>,
                        // Domain: x = U
                        flecsolve::multivariable_t<var::mass_density,
                          var::total_energy_density,
                          var::momentum_density>,
                        // Range: y = dU/dt
                        flecsolve::multivariable_t<var::mass_density,
                          var::total_energy_density,
                          var::momentum_density>> {
  using base = flecsolve::op::base<rhs_params<D>,
    flecsolve::multivariable_t<var::mass_density,
      var::total_energy_density,
      var::momentum_density>,
    flecsolve::multivariable_t<var::mass_density,
      var::total_energy_density,
      var::momentum_density>>;
  using base::params;

  explicit rhs_operator(rhs_params<D> p) : base(p) {}

  template<class Domain, class Range>
  void apply(Domain & x, Range & y) const {
    auto & s = *params.s;
    auto & sc = *params.sc;

    // Select the conserved variables from the input state x = U.
    auto & [rho, E, mom] = x;

    // Select the corresponding variables from the output RHS y = dU/dt.
    auto & [drho, dE, dmom] = y;

    /*
      Start from R = 0 before adding source and flux-divergence terms.
    */
    y.zero();

    // Primitive recovery
    annotation_rguard<tasks::hydro::conserved_to_primitive_task>(),
      sc.template execute<tasks::hydro::conserved_to_primitive<D>>(
        flecsi::exec::on,
        *s.m,
        rho.data.ref(),
        mom.data.ref(),
        E.data.ref(),
        s.prim.velocity(*s.m),
        s.prim.pressure(*s.m),
        s.prim.specific_internal_energy(*s.m),
        s.prim.sound_speed(*s.m),
        s.eos);

    // Apply boundary conditions.
    {
      annotation_rguard<tasks::apply_boundaries_task> guard;
      boundary_batch<D> bb(s.bcs);
      bb.add(s.cons.hydro.mass_density.name, rho.data.ref());
      bb.add(s.cons.hydro.total_energy_density.name, E.data.ref());
      bb.add(s.cons.hydro.momentum_density.name, mom.data.ref());
      bb.add(s.prim.pressure, *s.m);
      bb.add(s.prim.specific_internal_energy, *s.m);
      bb.add(s.prim.velocity, *s.m);
      bb.execute(sc, *s.m);
    }

    /*
      Explicit source contribution.

      This adds source terms to:
          d(rho*u)/dt
          d(rho*E)/dt
    */
    annotation_rguard<tasks::external_source_task>(),
      sc.template execute<tasks::external_source<D>>(flecsi::exec::on,
        *s.m,
        s.prim.velocity(*s.m),
        s.cons.hydro.mass_density(*s.m),
        s.src_t.hydro.gravity_accel(*s.m),
        dmom.data.ref(),
        dE.data.ref());

    for(std::size_t axis = 0; axis < D; axis++) {
      // Reconstruct primitive variables to faces.
      {
        auto scalar_v =
          std::vector{std::make_tuple(rho.data.ref(), s.f.hydro.r_face(s.m)),
            std::make_tuple(
              s.prim.specific_internal_energy(*s.m), s.f.hydro.e_face(s.m)),
            std::make_tuple(s.prim.sound_speed(*s.m), s.f.hydro.c_face(s.m)),
            std::make_tuple(s.prim.pressure(*s.m), s.f.hydro.p_face(s.m))};

        dispatch_reconstruct<D>(s.reconstructor,
          sc,
          axis,
          *s.m,
          scalar_v,
          std::vector{
            std::make_tuple(s.prim.velocity(*s.m), s.f.hydro.u_face(s.m))});
      }

      // Convert reconstructed primitive face states to conserved face states.
      annotation_rguard<tasks::hydro::get_conserved_task>(),
        sc.template execute<tasks::hydro::get_conserved<D>>(flecsi::exec::on,
          *s.m,
          s.f.hydro.r_face(s.m),
          s.f.hydro.u_face(s.m),
          s.f.hydro.e_face(s.m),
          s.f.hydro.ru_face(s.m),
          s.f.hydro.re_face(s.m));

      /*
        Compute numerical fluxes and add flux divergence to RHS:

            R += - div(F(U))

        In the current HARD implementation this task accumulates into the
        supplied RK/RHS storage.
      */
      annotation_rguard<tasks::hydro::compute_interface_fluxes_task>(),
        sc.template execute<tasks::hydro::compute_interface_fluxes<D>>(
          flecsi::exec::on,
          axis,
          *s.m,
          s.f.hydro.r_face(s.m),
          s.f.hydro.u_face(s.m),
          s.f.hydro.p_face(s.m),
          s.f.hydro.c_face(s.m),
          s.f.hydro.ru_face(s.m),
          s.f.hydro.re_face(s.m),
          s.rf.hydro.r_f(*s.m),
          s.rf.hydro.ru_f(*s.m),
          s.rf.hydro.re_f(*s.m),
          std::make_tuple(drho.data.ref(), dE.data.ref(), dmom.data.ref()),
          s.icst.gravity_const(*s.gt));
    }
  }
};

} // namespace hard::actions::flecsolve_hydro

#endif
