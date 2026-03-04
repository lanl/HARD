#ifndef HARD_MODULE_HYDRO_TASKS_INIT_HH
#define HARD_MODULE_HYDRO_TASKS_INIT_HH

#include "../constants.hh"
#include <cmath>
#include <flecsi/flog.hh>

namespace hard {

namespace tasks::init {

template<std::size_t D>
void inline initialize_gravity_force(flecsi::exec::cpu,
  typename field<vec<D>>::template accessor<wo, wo> gravity_force_a) noexcept {
  auto gf = gravity_force_a.span();
  // initialize gravity force with zero
  std::fill(gf.begin(), gf.end(), vec<D>{0.0});
} // initialize_gravity_force

template<std::size_t D>
void inline initialize_gravity_acc(
  typename single<vec<D>>::template accessor<wo> gravity_acc_a,
  vec<D> k) noexcept {
  (*gravity_acc_a) = k;
} // initialize_gravity_acc

void inline compute_dt_weighted(flecsi::exec::cpu,
  single<double>::accessor<ro> dt,
  single<double>::accessor<wo> dt_w,
  double tsg) noexcept {
  dt_w = *dt * tsg;
}

inline void
init_time(flecsi::exec::cpu,
  single<double>::accessor<wo> time,
  double vtime) noexcept {
  *time = vtime;
}

/*----------------------------------------------------------------------------*
  Temperature boundaries.
 *----------------------------------------------------------------------------*/

void inline set_t_boundary(flecsi::exec::cpu,
  field<double>::accessor<wo> t_boundary,
  std::vector<double> copy_values) noexcept {

  for(std::size_t i{0}; i < t_boundary.span().size(); i++) {
    t_boundary[i] = copy_values[i];
  }
} // t_boundary

/*----------------------------------------------------------------------------*
  Temperature unit conversion from eV (or other) to Kelvin
 *----------------------------------------------------------------------------*/

void inline convert_temperature(flecsi::exec::cpu,
  field<double>::accessor<rw> temperature,
  std::string const & unit) noexcept {

  assert((unit == "Kelvin" || unit == "eV") && "Unsupported temperature unit");

  double conversion_factor{};
  if(unit == "Kelvin") {
    return;
  }
  else if(unit == "eV") {
    conversion_factor = hard::constants::cgs::eV_to_K;
  }

  for(std::size_t i{0}; i < temperature.span().size(); i++) {
    temperature[i] *= conversion_factor;
  }
} // t_boundary

/*----------------------------------------------------------------------------*
  Set the average particle mass with a given mean molecular weight
 *----------------------------------------------------------------------------*/

void inline particle_mass(single<double>::accessor<wo> particle_mass_a,
  double mean_molecular_weight) noexcept {
  (*particle_mass_a) = mean_molecular_weight * constants::cgs::proton_mass;
} // particle_mass

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
  field<double>::accessor<wo, wo> // rEF_a
  ) noexcept {} // touch

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
  field<double>::accessor<wo, wo> // rEF_a
  ) noexcept {} // touch

} // namespace tasks::init
} // namespace hard

#endif // HARD_MODULE_HYDRO_TASKS_INIT_HH
