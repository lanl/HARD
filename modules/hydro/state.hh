#ifndef HARD_MODULE_HYDRO_STATE_HH
#define HARD_MODULE_HYDRO_STATE_HH

#include "../modules/spec/eos.hh"
#include "types.hh"

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

    static inline const typename single<
      typename mesh<D>::bmap>::template definition<global>
      bmap;
    static inline const single<double>::definition<global> particle_mass;
    static inline const field<double>::definition<global> time_boundary;
    static inline const field<double>::definition<global> temperature_boundary;

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
      static inline const field<double>::definition<mesh<D>, is::cells>
        mass_density;
      static inline const typename field<vec<D>>::template definition<mesh<D>,
        is::cells>
        momentum_density;
      static inline const field<double>::definition<mesh<D>, is::cells>
        total_energy_density;
    } hydro;

  } cons;

  // used in action, analyse, init
  struct primitives {
    // Primitives.
    static inline const typename field<vec<D>>::template definition<mesh<D>,
      is::cells>
      velocity; // u
    static inline const field<double>::definition<mesh<D>, is::cells>
      pressure; // p
    static inline const field<double>::definition<mesh<D>, is::cells>
      specific_internal_energy; // e
    static inline const field<double>::definition<mesh<D>, is::cells>
      sound_speed; // c
    static inline const field<double>::definition<mesh<D>, is::cells>
      temperature; // t
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

  // Storing dU/dt
  RK<D> rk_dt1;
  RK<D> rk_dt2;
  RK<D> rk_n;

}; // struct state

} // namespace hard::hydro

#endif // HARD_MODULE_HYDRO_STATE_HH
