#ifndef HARD_STATE_HH
#define HARD_STATE_HH

#include "spec/eos.hh"
#include "types.hh"
#ifdef USE_CATALYST
#include "catalyst/types.hh"
#endif // USE_CATALYST

#if defined(USE_FLECSOLVE) && USE_FLECSOLVE
#include "flecsolve/operators/core.hh"
#include "flecsolve/solvers/cg.hh"
#include "flecsolve/solvers/factory.hh"
#include "flecsolve/solvers/gmres.hh"
#include "flecsolve/vectors/topo_view.hh"
#endif

namespace hard {

/*----------------------------------------------------------------------------*
  Problem state.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct state {

  flecsi::future<double> dtmin_;

  /*--------------------------------------------------------------------------*
    Topology slots.
  *--------------------------------------------------------------------------*/

  eos::eos_wrapper eos;

  /*--------------------------------------------------------------------------*
    Topology slots.
   *--------------------------------------------------------------------------*/
  flecsi::topo::index::ptr ct; /* Color topology. */
  flecsi::topo::global::ptr gt;
  flecsi::topo::global::ptr dense_topology;

#ifdef USE_CATALYST
  flecsi::topo::index::ptr pt;
  /*----------------------------------------------------------------------------*
    Register the catalyst_data_structure on the "index" topology.
   *----------------------------------------------------------------------------*/
  static inline single<catalyst_attributes>::definition<flecsi::topo::index>
    catalyst_data;

#endif // USE_CATALYST

  // Deque for meshes in multigrid
  std::deque<typename mesh<D>::ptr> mh;

  // Define the fine grid
  typename mesh<D>::ptr & m = mh.emplace_back(typename mesh<D>::ptr());

  /*--------------------------------------------------------------------------*
    Global parameters.
   *--------------------------------------------------------------------------*/

  static inline const typename single<
    typename mesh<D>::bmap>::template definition<global>
    bmap;
#ifdef ENABLE_RADIATION
  static inline const single<double>::definition<global> kappa;
  static inline const single<std::size_t>::definition<global> limiter_id;
  static inline const single<std::size_t>::definition<global> closure_id;
#endif
  static inline const single<double>::definition<global> particle_mass;
  static inline const field<double>::definition<global> time_boundary;
  static inline const field<double>::definition<global> temperature_boundary;

  static inline const typename single<vec<D>>::template definition<global>
    gravity_acc;

  /*--------------------------------------------------------------------------*
    Color parameters (One per color using an index topology instance).
   *--------------------------------------------------------------------------*/

  /* Maximum characteristic speed for a color. */
  static inline const typename single<vec<D>>::template definition<index> lmax;

  static inline const single<double>::template definition<global> dt, t,
    dt_weighted;

  /*--------------------------------------------------------------------------*
    Mesh fields.
   *--------------------------------------------------------------------------*/

  // Conserved quantities.
  static inline const field<double>::definition<mesh<D>, is::cells>
    mass_density;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    momentum_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    total_energy_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    radiation_energy_density;

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

  // Faces.
  static inline const field<double>::definition<mesh<D>, is::cells> eTail;
  static inline const field<double>::definition<mesh<D>, is::cells> cTail;
  static inline const field<double>::definition<mesh<D>, is::cells> rTail;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    ruTail;
  static inline const field<double>::definition<mesh<D>, is::cells> rETail;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    uTail;
  static inline const field<double>::definition<mesh<D>, is::cells> pTail;
  static inline const field<double>::definition<mesh<D>, is::cells> EradTail;

  static inline const field<double>::definition<mesh<D>, is::cells> eHead;
  static inline const field<double>::definition<mesh<D>, is::cells> cHead;
  static inline const field<double>::definition<mesh<D>, is::cells> rHead;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    ruHead;
  static inline const field<double>::definition<mesh<D>, is::cells> rEHead;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    uHead;
  static inline const field<double>::definition<mesh<D>, is::cells> pHead;
  static inline const field<double>::definition<mesh<D>, is::cells> EradHead;

  // Riemann fluxes.
  static inline const field<double>::definition<mesh<D>, is::cells> rF;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    ruF;
  static inline const field<double>::definition<mesh<D>, is::cells> rEF;
  static inline const field<double>::definition<mesh<D>, is::cells> EradF;

  // Radiation pressure (P^{ij})
  static inline const typename field<spec::tensor<D,
    spec::tensor_rank::Two>>::template definition<mesh<D>, is::cells>
    radiation_pressure_tensor;

  // Variables related to the diffusion (multigrid) solver
  static inline dual_field<double, D> Esf; // Temp solution field in multigrid
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Uf; // Outer solution field
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Ef; // RHS of Au=f
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Ef_temp; // RHS of Au=f in multigrid precond
  static inline const typename field<stencil<D>>::template definition<mesh<D>,
    is::cells>
    Ew; // stencil weights
  static inline const field<double>::template definition<mesh<D>, is::cells>
    r; // rho
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Diff; // diffusion field.
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Df_x; // diffusion(Face) field.
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Df_y; // diffusion(Face) field.
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Df_z; // diffusion(Face) field.
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Resf; // Residual field.
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Errf; // Error field.

  // Number of levels that can be used.
  std::size_t max_num_levels{0};

  // Lowest level
  std::size_t lowest_level{0};

  // Highest level.
  std::size_t highest_level{0};

  // Pre smoothing
  std::size_t mg_pre{4};

  // Post smoothing
  std::size_t mg_post{4};

  // Cycles
  std::size_t mg_cycles{1};

  // Jacobi iterations in mg coarse grid
  std::size_t jacobi_iterations;

  // Gradient of a radiation energy density
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    gradient_rad_energy;

  // Magnitude of the gradient of the radiation energy density
  static inline const field<double>::definition<mesh<D>, is::cells>
    magnitude_gradient_rad_energy;

  // Gravity force
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    gravity_force;

  // Radiation force
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    radiation_force;

  // Dimensionless quantitiy, R
  static inline const field<double>::definition<mesh<D>, is::cells> R_value;

  // Flux limiter (standard/adaptive)
  static inline const field<double>::definition<mesh<D>, is::cells>
    lambda_bridge;

  // Eddington Factor
  static inline const field<double>::definition<mesh<D>, is::cells>
    eddington_factor;

  // Gradient of velocity
  static inline const typename field<spec::tensor<D,
    spec::tensor_rank::Two>>::template definition<mesh<D>, is::cells>
    velocity_gradient;

  // Storing dU/dt
  // In the future, it'd be worth consider grouping these into an one
  // std::tuple<> object per substep.
  //
  //  - For RK substep 1
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_mass_density;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    dt_momentum_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_total_energy_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_total_energy_density_implicit;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_radiation_energy_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_radiation_energy_density_implicit;
  //  - For RK substep 2
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_mass_density_2;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    dt_momentum_density_2;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_total_energy_density_2;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_total_energy_density_implicit_2;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_radiation_energy_density_2;
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_radiation_energy_density_implicit_2;

  // Storing variables at U^n during performing time update
  static inline const field<double>::definition<mesh<D>, is::cells>
    mass_density_n;
  static inline const typename field<vec<D>>::template definition<mesh<D>,
    is::cells>
    momentum_density_n;
  static inline const field<double>::definition<mesh<D>, is::cells>
    total_energy_density_n;
  static inline const field<double>::definition<mesh<D>, is::cells>
    radiation_energy_density_n;

  // FIXME Temporary variable for sections of the code only used for the
  // multigroup implementation (mostly boundaries)
  bool mg = false;

  /*--------------------------------------------------------------------------*
    FleCSolve
  *--------------------------------------------------------------------------*/

#if defined(USE_FLECSOLVE) && USE_FLECSOLVE
  flecsolve::bicgstab::settings solver_settings;
  bool flecsolve_coarse_grid;
  std::size_t nr_vcycles = 1;
#endif
  bool full_multigrid = false;

}; // struct state

} // namespace hard

#endif // HARD_STATE_HH
