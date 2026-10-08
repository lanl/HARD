#ifndef HARD_MODULE_HYDRO_STATE_HH
#define HARD_MODULE_HYDRO_STATE_HH

#include "../modules/spec/boundaries.hh"
#include "../modules/spec/eos.hh"
#include "types.hh"
#include <tuple>

namespace hard::hydro {

/*----------------------------------------------------------------------------*
  Problem state.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct state {

  /*--------------------------------------------------------------------------*
    EOS
  *--------------------------------------------------------------------------*/

  eos::eos_wrapper eos;

  /*--------------------------------------------------------------------------*
    Global parameters.
   *--------------------------------------------------------------------------*/

  struct initial_constants {

    static inline const single<double>::definition<global> particle_mass;

    static inline const typename single<vec<D>>::template definition<global>
      gravity_acc;

  } icst;

  /*--------------------------------------------------------------------------*
    Mesh fields.
   *--------------------------------------------------------------------------*/

  // used in action, analyse, init
  struct conserved {

    struct hydrodynamics {
      // Conserved quantities.
      static inline const tfield<double, mesh<D>> mass_density{"mass_density",
        tag::boundary};
      static inline const tfield<vec<D>, mesh<D>> momentum_density{
        "momentum_density",
        tag::boundary};
      static inline const tfield<double, mesh<D>> total_energy_density{
        "total_energy_density",
        tag::boundary};
    } hydro;

  } cons;

  // used in action, analyse, init
  struct primitives {
    // Primitives.
    static inline const tfield<vec<D>, mesh<D>> velocity{"velocity",
      tag::boundary}; // u
    static inline const tfield<double, mesh<D>> pressure{"pressure",
      tag::boundary}; // p
    static inline const tfield<double, mesh<D>> specific_internal_energy{
      "specific_internal_energy",
      tag::boundary}; // e
    static inline const tfield<double, mesh<D>> sound_speed{"sound_speed"}; // c
    static inline const tfield<double, mesh<D>> temperature{"temperature"}; // t
  } prim;

  // action
  struct faces_rl {
    struct hydro {
      faces<D> e_face;
      faces<D> c_face;
      faces<D> r_face;
      faces<D> re_face;
      faces<D> p_face;
      faces_vec<D> ru_face;
      faces_vec<D> u_face;
    } hydro;
  } f;

  // action
  struct rieman_fluxes {
    struct hydro {
      // Riemann fluxes.
      static inline const field<double>::definition<mesh<D>, is::cells> r_f;
      static inline const typename field<vec<D>>::template definition<mesh<D>,
        is::cells>
        ru_f;
      static inline const field<double>::definition<mesh<D>, is::cells> re_f;
    } hydro;
  } rf;

  // action, init
  struct source_term {
    struct hydrodynamics {
      // Gravity force
      static inline const typename field<vec<D>>::template definition<mesh<D>,
        is::cells>
        gravity_force;
    } hydro;
  } src_t;

  // Gradient of velocity
  static inline const typename field<spec::tensor<D,
    spec::tensor_rank::Two>>::template definition<mesh<D>, is::cells>
    velocity_gradient;

  /*--------------------------------------------------------------------------*
    Tagged fields and boundary conditions.
    *--------------------------------------------------------------------------*/

  auto fields() const {
    return std::tie(cons.hydro.mass_density,
      cons.hydro.momentum_density,
      cons.hydro.total_energy_density,
      prim.velocity,
      prim.pressure,
      prim.specific_internal_energy,
      prim.sound_speed,
      prim.temperature);
  }

  spec::boundaries<D> bcs;

  // Storing dU/dt
  // FIXME: We will not need this anymore
  RK<D> rk_dt1;
  RK<D> rk_dt2;
  RK<D> rk_n;

  spec::rc::reconstructor reconstructor = spec::rc::weno5z;
  spec::slv::solver use_solver = spec::slv::heun;

}; // struct state

} // namespace hard::hydro

#endif // HARD_MODULE_HYDRO_STATE_HH
