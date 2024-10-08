#ifndef HARD_TASKS_INIT_HH
#define HARD_TASKS_INIT_HH

#include "../constants.hh"
#include "../types.hh"
#include <cmath>
#include <flecsi/flog.hh>

namespace hard {

namespace tasks::init {

void inline compute_dt_weighted(single<double>::accessor<ro> dt,
  single<double>::accessor<wo> dt_w,
  double tsg) {
  dt_w = *dt * tsg;
}

inline void
init_time(single<double>::accessor<wo> time, double vtime) {
  *time = vtime;
}

/*----------------------------------------------------------------------------*
  Ideal gas parameter.
 *----------------------------------------------------------------------------*/

void inline gamma(single<double>::accessor<wo> gamma_a, double g) {
  (*gamma_a) = g;
} // gamma

/*----------------------------------------------------------------------------*
  Opacity parameter.
 *----------------------------------------------------------------------------*/

void inline kappa(single<double>::accessor<wo> kappa_a, double k) {
  (*kappa_a) = k;
} // kappa

/*----------------------------------------------------------------------------*
  Set the average particle mass with a given mean molecular weight
 *----------------------------------------------------------------------------*/

void inline particle_mass(single<double>::accessor<wo> particle_mass_a,
  double mean_molecular_weight) {
  (*particle_mass_a) = mean_molecular_weight * constants::cgs::proton_mass;
} // kappa

/*----------------------------------------------------------------------------*
  Fake initialization tasks to avoid legion errors.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
inline void
touch(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, wo> r_a,
  typename field<vec<D>>::template accessor<wo, wo> ru_a,
  field<double>::accessor<wo, wo> rE_a,

  typename field<vec<D>>::template accessor<wo, wo> u_a,
  field<double>::accessor<wo, wo> p_a,
  field<double>::accessor<wo, wo> rTail_a,
  typename field<vec<D>>::template accessor<wo, wo> ruTail_a,
  field<double>::accessor<wo, wo> rETail_a,
  typename field<vec<D>>::template accessor<wo, wo> uTail_a,
  field<double>::accessor<wo, wo> pTail_a,

  field<double>::accessor<wo, wo> rHead_a,
  typename field<vec<D>>::template accessor<wo, wo> ruHead_a,
  field<double>::accessor<wo, wo> rEHead_a,
  typename field<vec<D>>::template accessor<wo, wo> uHead_a,
  field<double>::accessor<wo, wo> pHead_a,

  field<double>::accessor<wo, wo> rF_a,
  typename field<vec<D>>::template accessor<wo, wo> ruF_a,
  field<double>::accessor<wo, wo> rEF_a
#ifndef DISABLE_RADIATION
  ,
  field<double>::accessor<wo, wo> Erad_a,
  field<double>::accessor<wo, wo> EradTail_a,
  field<double>::accessor<wo, wo> EradHead_a,
  field<double>::accessor<wo, wo> EradF_a,
  field<double>::accessor<wo, wo> Esf_a,
  field<double>::accessor<wo, wo> Ef_a,
  typename field<stencil<D>>::template accessor<wo, wo> Ew_a,
  field<double>::accessor<wo, wo> Df_a,
  field<double>::accessor<wo, wo> Resf_a,
  field<double>::accessor<wo, wo> Errf_a,
  //
  typename field<vec<D>>::template accessor<wo, wo> gradient_rad_energy_a,
  field<double>::accessor<wo, wo> magnitude_gradient_rad_energy_a,
  typename field<vec<D>>::template accessor<wo, wo> radiation_force_a,
  field<double>::accessor<wo, wo> R_value_a,
  field<double>::accessor<wo, wo> lambda_bridge_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    wo> velocity_gradient_a
#endif
) {
} // touch

template<std::size_t D>
inline void
touch1(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, wo> r_a,
  typename field<vec<D>>::template accessor<wo, wo> ru_a,
  field<double>::accessor<wo, wo> rE_a,

  typename field<vec<D>>::template accessor<wo, wo> u_a,
  field<double>::accessor<wo, wo> p_a,
  field<double>::accessor<wo, wo> rTail_a,
  typename field<vec<D>>::template accessor<wo, wo> ruTail_a,
  field<double>::accessor<wo, wo> rETail_a,
  typename field<vec<D>>::template accessor<wo, wo> uTail_a,
  field<double>::accessor<wo, wo> pTail_a,

  field<double>::accessor<wo, wo> rHead_a,
  typename field<vec<D>>::template accessor<wo, wo> ruHead_a,
  field<double>::accessor<wo, wo> rEHead_a,
  typename field<vec<D>>::template accessor<wo, wo> uHead_a,
  field<double>::accessor<wo, wo> pHead_a,

  field<double>::accessor<wo, wo> rF_a,
  typename field<vec<D>>::template accessor<wo, wo> ruF_a,
  field<double>::accessor<wo, wo> rEF_a
#ifndef DISABLE_RADIATION
  ,
  field<double>::accessor<wo, wo> Erad_a,
  field<double>::accessor<wo, wo> EradTail_a,
  field<double>::accessor<wo, wo> EradHead_a,
  field<double>::accessor<wo, wo> EradF_a,
  field<double>::accessor<wo, wo> Esf_a,
  field<double>::accessor<wo, wo> Ef_a,
  typename field<stencil<D>>::template accessor<wo, wo> Ew_a,
  field<double>::accessor<wo, wo> Df_a,
  field<double>::accessor<wo, wo> Resf_a,
  field<double>::accessor<wo, wo> Errf_a,
  //
  typename field<vec<D>>::template accessor<wo, wo> gradient_rad_energy_a,
  field<double>::accessor<wo, wo> magnitude_gradient_rad_energy_a,
  typename field<vec<D>>::template accessor<wo, wo> radiation_force_a,
  field<double>::accessor<wo, wo> R_value_a,
  field<double>::accessor<wo, wo> lambda_bridge_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    wo> velocity_gradient_a
#endif
) {
} // touch

/*----------------------------------------------------------------------------*
  Boundary input translation.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
mesh<D>::periodic_axes
boundaries(
  typename single<typename mesh<D>::bmap>::template accessor<wo> bmap_a,
  std::array<std::array<bd::boundary_type, 2>, D> bnds) {
  auto & bmap = *bmap_a;
  constexpr auto periodic = mesh<D>::boundary_type::periodic;

  typename mesh<D>::periodic_axes p;
  p[ax::x] = false;
  if constexpr(D == 2 || D == 3) {
    p[ax::y] = false;
  }
  else if constexpr(D == 3) {
    p[ax::z] = false;
  } // if

  bmap[ax::x][bd::low] = bnds[ax::x][bd::low];
  bmap[ax::x][bd::high] = bnds[ax::x][bd::high];
  p[ax::x] =
    bnds[ax::x][bd::low] == periodic && bnds[ax::x][bd::high] == periodic;

  if constexpr(D == 2 || D == 3) {
    bmap[ax::y][bd::low] = bnds[ax::y][bd::low];
    bmap[ax::y][bd::high] = bnds[ax::y][bd::high];
    p[ax::y] =
      bnds[ax::y][bd::low] == periodic && bnds[ax::y][bd::high] == periodic;
  } // if

  if constexpr(D == 3) {
    bmap[ax::z][bd::low] = bnds[ax::z][bd::low];
    bmap[ax::z][bd::high] = bnds[ax::z][bd::high];
    p[ax::z] =
      bnds[ax::z][bd::low] == periodic && bnds[ax::z][bd::high] == periodic;
  } // if

  return p;
} // boundaries

} // namespace tasks::init
} // namespace hard

#endif // HARD_TASKS_INIT_HH
