#ifndef HARD_MODULE_HYDRO_TIME_DERIVATIVE_HH
#define HARD_MODULE_HYDRO_TIME_DERIVATIVE_HH

#include "../modules/hydro/numerical_algorithms/time_stepper.hh"
#include <cstddef>
#include <flecsi/utilities.hh>

namespace hard::tasks::hydro {

template<std::size_t D>
double
update_dtmin(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  flecsi::future<double> lmax_f) noexcept {

  double lmax = lmax_f.get();
  if constexpr(D == 1) {
    return m.template delta<ax::x>() / lmax;
  }
  else if constexpr(D == 2) {
    return std::min(m.template delta<ax::x>(), m.template delta<ax::y>()) /
           (D * lmax);
  }
  else {
    return std::min(m.template delta<ax::x>(),
             std::min(m.template delta<ax::y>(), m.template delta<ax::z>())) /
           (D * lmax);
  } // if
} // update_dtmin

template<std::size_t Dim>
void
set_dudt_to_zero(flecsi::exec::accelerator s,
  std::vector<field<double>::accessor<wo, na>> rk_dt_v_a,
  std::vector<typename field<vec<Dim>>::accessor<wo, na>>
    rk_dt_vec_v_a) noexcept {

  for(auto rk_dt_a : rk_dt_v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, rk_dt_a.span().size())) {
      rk_dt_a(i) = 0.0;
    }; // forall
  }

  for(auto rk_dt_a : rk_dt_vec_v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, rk_dt_a.span().size())) {
      rk_dt_a(i) = vec<Dim>(0.0);
    }; // forall
  }
}

//
// Store the evolved variables U^n at t=t^n into temporary space before running
// RK substeps.
//
template<std::size_t Dim>
void
store_current_state(flecsi::exec::accelerator s,
  std::vector<field<double>::accessor<ro, na>> from_v_a,
  std::vector<typename field<vec<Dim>>::accessor<ro, na>> from_vec_v_a,
  std::vector<field<double>::accessor<wo, na>> to_v_a,
  std::vector<typename field<vec<Dim>>::accessor<wo, na>> to_vec_v_a) noexcept {

  for(std::size_t i = 0; i < from_v_a.size(); ++i) {
    auto && from_a = from_v_a[i];
    auto && to_a = to_v_a[i];
    s.executor().forall(i, flecsi::util::iota_view({}, from_a.span().size())) {
      to_a(i) = from_a(i);
    }; // forall
  }

  for(std::size_t i = 0; i < from_vec_v_a.size(); ++i) {
    auto && from_a = from_vec_v_a[i];
    auto && to_a = to_vec_v_a[i];
    s.executor().forall(i, flecsi::util::iota_view({}, from_a.span().size())) {
      to_a(i) = from_a(i);
    }; // forall
  }
}

//
// Used for the RK substeps
//
template<std::size_t Dim>
void
update_u(flecsi::exec::accelerator s,
  single<double>::accessor<ro> dt_a,
  // U^n we want to update
  std::vector<field<double>::accessor<rw, na>> rk_n_v_a,
  std::vector<typename field<vec<Dim>>::accessor<rw, na>> rk_n_vec_v_a,
  // Time derivatives for the state U^1
  std::vector<field<double>::accessor<ro, na>> rk_dt_v_a,
  std::vector<typename field<vec<Dim>>::accessor<ro, na>>
    rk_dt_vec_v_a) noexcept {

  auto h = *dt_a;

  // Scalar
  for(std::size_t i = 0; i < rk_n_v_a.size(); ++i) {
    auto && rk_n_a = rk_n_v_a[i];
    auto && rk_dt_a = rk_dt_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_n_a.span().size())) {
      rk_n_a(j) += h * rk_dt_a(j);
    }; // forall
  }

  // Vec<D>
  for(std::size_t i = 0; i < rk_n_vec_v_a.size(); ++i) {
    auto && rk_n_a = rk_n_vec_v_a[i];
    auto && rk_dt_a = rk_dt_vec_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_n_a.span().size())) {
      rk_n_a(j) += h * rk_dt_a(j);
    }; // forall
  }
}

template<std::size_t Dim>
void
update_u_scalar(flecsi::exec::accelerator s,
  single<double>::accessor<ro> dt_a,
  // U^n we want to update
  std::vector<field<double>::accessor<rw, na>> rk_n_v_a,
  // Time derivatives for the state U^1
  std::vector<field<double>::accessor<ro, na>> rk_dt_v_a) noexcept {

  auto h = *dt_a;

  // Scalar
  for(std::size_t i = 0; i < rk_n_v_a.size(); ++i) {
    auto && rk_n_a = rk_n_v_a[i];
    auto && rk_dt_a = rk_dt_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_n_a.span().size())) {
      rk_n_a(j) += h * rk_dt_a(j);
    }; // forall
  }
}

//
// Used for the RK stage updates
//
template<std::size_t Dim>
void
update_u_stage(flecsi::exec::cpu s,
  single<double>::accessor<ro> dt_a,
  // U^n we want to update
  // U^n we want to update
  std::vector<field<double>::accessor<ro, na>> rk_n_v_a,
  std::vector<typename field<vec<Dim>>::accessor<ro, na>> rk_n_vec_v_a,
  // Time derivatives for the state U^1
  std::vector<field<double>::accessor<ro, na>> rk_dt_v_a,
  std::vector<typename field<vec<Dim>>::accessor<ro, na>> rk_dt_vec_v_a,
  // U^n+1 updated after stage
  std::vector<field<double>::accessor<rw, na>> rk_new_v_a,
  std::vector<typename field<vec<Dim>>::accessor<rw, na>>
    rk_new_vec_v_a) noexcept {

  auto h = *dt_a;
  // Scalar
  for(std::size_t i = 0; i < rk_n_v_a.size(); ++i) {
    auto && rk_n_a = rk_n_v_a[i];
    auto && rk_dt_a = rk_dt_v_a[i];
    auto && rk_new_a = rk_new_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_n_a.span().size())) {
      rk_new_a(j) += rk_n_a(i) + h * rk_dt_a(j);
    }; // forall
  }

  // Vec<D>
  for(std::size_t i = 0; i < rk_n_vec_v_a.size(); ++i) {
    auto && rk_n_a = rk_n_vec_v_a[i];
    auto && rk_dt_a = rk_dt_vec_v_a[i];
    auto && rk_new_a = rk_new_vec_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_n_a.span().size())) {
      rk_new_a(j) += rk_n_a(j) + h * rk_dt_a(j);
    }; // forall
  }
}

template<std::size_t Dim>
void
add_k1_k2(flecsi::exec::accelerator s,
  // RK1
  std::vector<field<double>::accessor<rw, na>> rk_1_v_a,
  std::vector<typename field<vec<Dim>>::accessor<rw, na>> rk_1_vec_v_a,
  // RK2
  std::vector<field<double>::accessor<ro, na>> rk_2_v_a,
  std::vector<typename field<vec<Dim>>::accessor<ro, na>>
    rk_2_vec_v_a) noexcept {

  // Scalar
  for(std::size_t i = 0; i < rk_1_v_a.size(); ++i) {
    auto && rk_1_a = rk_1_v_a[i];
    auto && rk_2_a = rk_2_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_1_a.span().size())) {
      rk_1_a(j) = (rk_1_a(j) + rk_2_a(j)) * 0.5;
    }; // forall
  }

  // Vec<D>
  for(std::size_t i = 0; i < rk_1_vec_v_a.size(); ++i) {
    auto && rk_1_a = rk_1_vec_v_a[i];
    auto && rk_2_a = rk_2_vec_v_a[i];
    s.executor().forall(j, flecsi::util::iota_view({}, rk_1_a.span().size())) {
      rk_1_a(j) = (rk_1_a(j) + rk_2_a(j)) * 0.5;
    }; // forall
  }
}

} // namespace hard::tasks::hydro

#endif // HARD_MODULE_HYDRO_TIME_DERIVATIVE_HH
