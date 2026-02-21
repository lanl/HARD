#ifndef HARD_MODULE_RAD_STATE_HH
#define HARD_MODULE_RAD_STATE_HH

#include "../modules/spec/eos.hh"

#include "flecsolve/operators/core.hh"
#include "flecsolve/solvers/cg.hh"
#include "flecsolve/solvers/factory.hh"
#include "flecsolve/solvers/gmres.hh"
#include "flecsolve/vectors/topo_view.hh"

namespace hard::rad {

/*----------------------------------------------------------------------------*
  Problem state.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct state {

  struct rad {

    /*--------------------------------------------------------------------------*
      Global parameters.
      *--------------------------------------------------------------------------*/
    struct initial_constants {
      static inline const single<double>::definition<global> kappa;
      static inline const single<std::size_t>::definition<global> limiter_id;
      static inline const single<std::size_t>::definition<global> closure_id;
    } icst;

    /*--------------------------------------------------------------------------*
      Mesh fields.
      *--------------------------------------------------------------------------*/

    struct conservatives {
      static inline const field<double>::definition<mesh<D>, is::cells>
        radiation_energy_density;
    } cons;

    struct faces_rl {
      faces<D> EradFace;
    } f;

    struct rieman_fluxes {
      static inline const field<double>::definition<mesh<D>, is::cells> EradF;
    } rf;

    struct multigrid {
      // Variables related to the diffusion (multigrid) solver
      static inline dual_field<double, D>
        Esf; // Temp solution field in multigrid
      static inline const field<double>::template definition<mesh<D>, is::cells>
        Uf; // Outer solution field
      static inline const field<double>::template definition<mesh<D>, is::cells>
        Ef; // RHS of Au=f
      static inline const field<double>::template definition<mesh<D>, is::cells>
        Ef_temp; // RHS of Au=f in multigrid precond
      static inline const typename field<stencil<D>>::template definition<
        mesh<D>,
        is::cells>
        Ew; // stencil weights
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

      // Pre smoothing
      std::size_t mg_pre{4};

      // Post smoothing
      std::size_t mg_post{4};

      // Cycles
      std::size_t mg_cycles{1};

      // Jacobi iterations in mg coarse grid
      std::size_t jacobi_iterations;

      /*--------------------------------------------------------------------------*
        FleCSolve
        *--------------------------------------------------------------------------*/

      flecsolve::bicgstab::settings solver_settings;
      bool flecsolve_coarse_grid;
      std::size_t nr_vcycles = 1;
      bool full_multigrid = false;

    } mgr;

    struct source_term {
      // Radiation force
      static inline const typename field<vec<D>>::template definition<mesh<D>,
        is::cells>
        radiation_force;
      // Radiation pressure (P^{ij})
      static inline const typename field<spec::tensor<D,
        spec::tensor_rank::Two>>::template definition<mesh<D>, is::cells>
        radiation_pressure_tensor;
    } src_t;

    struct limiter {

      // Gradient of a radiation energy density
      static inline const typename field<vec<D>>::template definition<mesh<D>,
        is::cells>
        gradient_rad_energy;

      // Magnitude of the gradient of the radiation energy density
      static inline const field<double>::definition<mesh<D>, is::cells>
        magnitude_gradient_rad_energy;

      // Dimensionless quantitiy, R
      static inline const field<double>::definition<mesh<D>, is::cells> R_value;

      // Flux limiter (standard/adaptive)
      static inline const field<double>::definition<mesh<D>, is::cells>
        lambda_bridge;

      // Eddington Factor
      static inline const field<double>::definition<mesh<D>, is::cells>
        eddington_factor;

    } limiter;

    field<double>::definition<mesh<D>, is::cells> dt_radiation_energy_density_1;
    field<double>::definition<mesh<D>, is::cells> dt_radiation_energy_density_2;
    field<double>::definition<mesh<D>, is::cells> dt_radiation_energy_density_n;

  } rad; // struct rad

}; // struct state

} // namespace hard::rad

#endif // HARD_MODULE_RAD_STATE_HH
