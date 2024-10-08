#ifndef HARD_STATE_HH
#define HARD_STATE_HH

#include "types.hh"
#ifdef USE_CATALYST
#include "catalyst/types.hh"
#endif // USE_CATALYST

namespace hard {

/*----------------------------------------------------------------------------*
  Global parameters.
 *----------------------------------------------------------------------------*/

inline const single<double>::definition<flecsi::topo::global> gamma;
inline const single<double>::definition<flecsi::topo::global> kappa;

inline const single<double>::definition<flecsi::topo::global> particle_mass;

/*----------------------------------------------------------------------------*
  Problem state.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct state {

  flecsi::future<double> dtmin_;

  /*--------------------------------------------------------------------------*
    Topology slots.
   *--------------------------------------------------------------------------*/
  index::slot ct; /* Color topology. */
  flecsi::topo::global::slot gt;

  // Vector for meshes in multigrid
  std::vector<std::unique_ptr<typename mesh<D>::slot>> mh;

  // Define the fine grid
  typename mesh<D>::slot & m =
    *mh.emplace_back(std::make_unique<typename mesh<D>::slot>());

  /*--------------------------------------------------------------------------*
    Global parameters.
   *--------------------------------------------------------------------------*/

  static inline const single<typename mesh<D>::bmap>::template definition<
    flecsi::topo::global>
    bmap;

  /*--------------------------------------------------------------------------*
    Color parameters (One per color using an index topology instance).
   *--------------------------------------------------------------------------*/

  /* Maximum characteristic speed for a color. */
  static inline const single<vec<D>>::template definition<flecsi::topo::index>
    lmax;

  static inline const single<double>::template definition<flecsi::topo::global>
    dt, t, dt_weighted;

  /*--------------------------------------------------------------------------*
    Mesh fields.
   *--------------------------------------------------------------------------*/

  // Conserved quantities.
  static inline const field<double>::definition<mesh<D>, is::cells>
    mass_density;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    momentum_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    total_energy_density;
  static inline const field<double>::definition<mesh<D>, is::cells>
    radiation_energy_density;

  // Primitives.
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    velocity;
  static inline const field<double>::definition<mesh<D>, is::cells> pressure;

  // Faces.
  static inline const field<double>::definition<mesh<D>, is::cells> rTail;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    ruTail;
  static inline const field<double>::definition<mesh<D>, is::cells> rETail;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    uTail;
  static inline const field<double>::definition<mesh<D>, is::cells> pTail;
  static inline const field<double>::definition<mesh<D>, is::cells> EradTail;

  static inline const field<double>::definition<mesh<D>, is::cells> rHead;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    ruHead;
  static inline const field<double>::definition<mesh<D>, is::cells> rEHead;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    uHead;
  static inline const field<double>::definition<mesh<D>, is::cells> pHead;
  static inline const field<double>::definition<mesh<D>, is::cells> EradHead;

  // Riemann fluxes.
  static inline const field<double>::definition<mesh<D>, is::cells> rF;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    ruF;
  static inline const field<double>::definition<mesh<D>, is::cells> rEF;
  static inline const field<double>::definition<mesh<D>, is::cells> EradF;

  // Radiation pressure (P^{ij})
  static inline const field<spec::tensor<D,
    spec::tensor_rank::Two>>::template definition<mesh<D>, is::cells>
    radiation_pressure_tensor;

  // Variables related to the diffusion (multigrid) solver
  static inline dual_field<double, D> Esf; // solution field
  static inline const field<double>::template definition<mesh<D>, is::cells>
    Ef; // RHS of Au=f
  static inline const field<stencil<D>>::template definition<mesh<D>, is::cells>
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

  // Gradient of a radiation energy density
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    gradient_rad_energy;

  // Magnitude of the gradient of the radiation energy density
  static inline const field<double>::definition<mesh<D>, is::cells>
    magnitude_gradient_rad_energy;

  // Radiation force
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    radiation_force;

  // Dimensionless quantitiy, R
  static inline const field<double>::definition<mesh<D>, is::cells> R_value;

  // Flux limiter
  static inline const field<double>::definition<mesh<D>, is::cells>
    lambda_bridge;

  // Gradient of velocity
  static inline const field<spec::tensor<D,
    spec::tensor_rank::Two>>::template definition<mesh<D>, is::cells>
    velocity_gradient;

  // Storing dU/dt
  // In the future, it'd be worth consider grouping these into an one
  // std::tuple<> object per substep.
  //
  //  - For RK substep 1
  static inline const field<double>::definition<mesh<D>, is::cells>
    dt_mass_density;
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
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
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
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
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    momentum_density_n;
  static inline const field<double>::definition<mesh<D>, is::cells>
    total_energy_density_n;
  static inline const field<double>::definition<mesh<D>, is::cells>
    radiation_energy_density_n;

}; // struct state

#ifdef USE_CATALYST

/*----------------------------------------------------------------------------*
  Convenience reference of the flecsi pre-defined process_topology slot. The
  process_topology slot has the same number of colors as there are processes.
  This means that it can be used as an "mpi" topology, i.e., the number of
  colors will always match the number of mpi ranks. In this example each
  process/color will have its own instance of a type that is registered on it.
 *----------------------------------------------------------------------------*/
inline flecsi::topo::index::slot & pt = flecsi::process_topology;

/*----------------------------------------------------------------------------*
  Register the catalyst_data_structure on the "index" topology.
 *----------------------------------------------------------------------------*/
inline const single<catalyst_attributes>::definition<flecsi::topo::index>
  catalyst_data;

#endif // USE_CATALYST

} // namespace hard

#endif // HARD_STATE_HH
