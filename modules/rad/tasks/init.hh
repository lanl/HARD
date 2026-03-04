#ifndef HARD_MODULE_RAD_TASKS_INIT_HH
#define HARD_MODULE_RAD_TASKS_INIT_HH

#include "../constants.hh"
#include "types.hh"
#include <cmath>
#include <flecsi/flog.hh>

namespace hard {

namespace tasks::init {

/*----------------------------------------------------------------------------*
  Adaptive Limiter Check, Closure ID, Limiter ID
 *----------------------------------------------------------------------------*/

// Simple opacity that is used across all
void inline closure_id(single<std::size_t>::accessor<wo> closure_id_a,
  std::size_t clid) noexcept {
  (*closure_id_a) = clid;
} // closure_id

// Simple opacity that is used across all
void inline limiter_id(single<std::size_t>::accessor<wo> limiter_id_a,
  std::size_t lmid) noexcept {
  (*limiter_id_a) = lmid;
} // limiter_id

/*----------------------------------------------------------------------------*
  Opacity parameter.
 *----------------------------------------------------------------------------*/

// Simple opacity that is used across all
void inline kappa(single<double>::accessor<wo> kappa_a, double k) {
  (*kappa_a) = k;
} // kappa

/*----------------------------------------------------------------------------*
  Fake initialization tasks to avoid legion errors.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
inline void
touch(typename mesh<D>::template accessor<ro>, // m,
  field<double>::accessor<wo, wo>, // r_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ru_a,
  field<double>::accessor<wo, wo>, // rE_a,

  typename field<vec<D>>::template accessor<wo, wo>, // u_a,
  field<double>::accessor<wo, wo>, // p_a,
  field<double>::accessor<wo, wo>, // rLeft_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ruLeft_a,
  field<double>::accessor<wo, wo>, // rELeft_a,
  typename field<vec<D>>::template accessor<wo, wo>, // uLeft_a,
  field<double>::accessor<wo, wo>, // pLeft_a,

  field<double>::accessor<wo, wo>, // rRight_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ruRight_a,
  field<double>::accessor<wo, wo>, // rERight_a,
  typename field<vec<D>>::template accessor<wo, wo>, // uRight_a,
  field<double>::accessor<wo, wo>, // pRight_a,

  field<double>::accessor<wo, wo>, // rF_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ruF_a,
  field<double>::accessor<wo, wo>, // rEF_a
  field<double>::accessor<wo, wo>, // Erad_a,
  field<double>::accessor<wo, wo>, // EradLeft_a,
  field<double>::accessor<wo, wo>, // EradRight_a,
  field<double>::accessor<wo, wo>, // EradF_a,
  field<double>::accessor<wo, wo>, // Esf_a,
  field<double>::accessor<wo, wo>, // Ef_a,
  typename field<stencil<D>>::template accessor<wo, wo>, // Ew_a,
  field<double>::accessor<wo, wo>, // Df_a,
  field<double>::accessor<wo, wo>, // Resf_a,
  field<double>::accessor<wo, wo>, // Errf_a,
  //
  typename field<vec<D>>::template accessor<wo, wo>, // gradient_rad_energy_a,
  field<double>::accessor<wo, wo>, // magnitude_gradient_rad_energy_a,
  typename field<vec<D>>::template accessor<wo, wo>, // radiation_force_a,
  field<double>::accessor<wo, wo>, // R_value_a,
  field<double>::accessor<wo, wo>, // lambda_bridge_a,
  field<double>::accessor<wo, wo>, // eddington_factor_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    wo> // velocity_gradient_a
) {} // touch

template<std::size_t D>
inline void
touch1(typename mesh<D>::template accessor<ro>, // m,
  field<double>::accessor<wo, wo>, // r_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ru_a,
  field<double>::accessor<wo, wo>, // rE_a,

  typename field<vec<D>>::template accessor<wo, wo>, // u_a,
  field<double>::accessor<wo, wo>, // p_a,
  field<double>::accessor<wo, wo>, // rLeft_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ruLeft_a,
  field<double>::accessor<wo, wo>, // rELeft_a,
  typename field<vec<D>>::template accessor<wo, wo>, // uLeft_a,
  field<double>::accessor<wo, wo>, // pLeft_a,

  field<double>::accessor<wo, wo>, // rRight_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ruRight_a,
  field<double>::accessor<wo, wo>, // rERight_a,
  typename field<vec<D>>::template accessor<wo, wo>, // uRight_a,
  field<double>::accessor<wo, wo>, // pRight_a,

  field<double>::accessor<wo, wo>, // rF_a,
  typename field<vec<D>>::template accessor<wo, wo>, // ruF_a,
  field<double>::accessor<wo, wo>, // rEF_a
  field<double>::accessor<wo, wo>, // Erad_a,
  field<double>::accessor<wo, wo>, // EradLeft_a,
  field<double>::accessor<wo, wo>, // EradRight_a,
  field<double>::accessor<wo, wo>, // EradF_a,
  field<double>::accessor<wo, wo>, // Esf_a,
  field<double>::accessor<wo, wo>, // Ef_a,
  typename field<stencil<D>>::template accessor<wo, wo>, // Ew_a,
  field<double>::accessor<wo, wo>, // Df_a,
  field<double>::accessor<wo, wo>, // Resf_a,
  field<double>::accessor<wo, wo>, // Errf_a,
  //
  typename field<vec<D>>::template accessor<wo, wo>, // gradient_rad_energy_a,
  field<double>::accessor<wo, wo>, // magnitude_gradient_rad_energy_a,
  typename field<vec<D>>::template accessor<wo, wo>, // radiation_force_a,
  field<double>::accessor<wo, wo>, // R_value_a,
  field<double>::accessor<wo, wo>, // lambda_bridge_a,
  field<double>::accessor<wo, wo>, // eddington_factor_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    wo> // velocity_gradient_a

) {} // touch

} // namespace tasks::init
} // namespace hard

#endif // HARD_MODULE_RAD_TASKS_INIT_HH
