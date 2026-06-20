#ifndef HARD_MODULE_HYDRO_RHS_HH
#define HARD_MODULE_HYDRO_RHS_HH

#include "../modules/hydro/numerical_algorithms/time_stepper.hh"
#include <flecsi/utilities.hh>

namespace hard::tasks::hydro {

struct update_dtmin_task : annotation_task_region {
  inline static const std::string name{"update_dtmin"};
};

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

struct set_dudt_to_zero_task : annotation_task_region {
  inline static const std::string name{"set_dudt_to_zero"};
};

template<std::size_t Dim>
void
set_dudt_to_zero(flecsi::exec::accelerator s,
  std::vector<field<double>::accessor<wo, na>> rk_dt_v_a,
  std::vector<typename vector_field_double<Dim>::template accessor<wo, na>>
    rk_dt_vec_v_a) noexcept {

  for(auto rk_dt_a : rk_dt_v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, rk_dt_a.span().size())) {
      rk_dt_a(i) = 0.0;
    }; // forall
  }

  for(auto rk_dt_a : rk_dt_vec_v_a) {
    auto gpu_rk_dt_a = vector_field_double<Dim>::to_accelerator(rk_dt_a);
    for(int axis = 0; axis < rk_dt_a.size(); axis++) {
      s.executor().forall(
        i, flecsi::util::iota_view({}, gpu_rk_dt_a[axis].span().size())) {
        gpu_rk_dt_a[axis](i) = 0.0;
      }; // forall
    }
  }
}

struct update_u_task : annotation_task_region {
  inline static const std::string name{"update_u"};
};

//
// Used for the RK substeps
//
template<std::size_t Dim>
void
update_u(flecsi::exec::accelerator s,
  single<double>::accessor<ro> h_a,
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    field<double>::accessor<rw, na>>> v_a,
  std::vector<
    std::tuple<typename vector_field_double<Dim>::template accessor<ro, na>,
      typename vector_field_double<Dim>::template accessor<rw, na>>>
    vec_v_a) noexcept {

  // Scalar
  for(auto & [dt_a, a] : v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, a.span().size())) {
      a(i) += *h_a * dt_a(i);
    }; // forall
  }

  // Vec<D>
  for(auto & [vec_dt_a, vec_a] : vec_v_a) {
    auto gpu_vec_dt_a = vector_field_double<Dim>::to_accelerator(vec_dt_a);
    auto gpu_vec_a = vector_field_double<Dim>::to_accelerator(vec_a);
    for(int axis = 0; axis < gpu_vec_a.size(); axis++) {
      s.executor().forall(
        i, flecsi::util::iota_view({}, gpu_vec_dt_a[axis].span().size())) {
        gpu_vec_a[axis](i) += *h_a * gpu_vec_dt_a[axis](i);
      }; // forall
    }
  }
}

struct update_u_scalar_task : annotation_task_region {
  inline static const std::string name{"update_u_scalar"};
};

template<std::size_t Dim>
void
update_u_scalar(flecsi::exec::accelerator s,
  single<double>::accessor<ro> h_a,
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    field<double>::accessor<rw, na>>> v_a) noexcept {

  // Scalar
  for(auto & [dt_a, a] : v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, a.span().size())) {
      a(i) += *h_a * dt_a(i);
    }; // forall
  }
}

struct update_u_stage_task : annotation_task_region {
  inline static const std::string name{"update_u_stage"};
};

//
// Used for the RK stage updates
//
template<std::size_t Dim>
void
update_u_stage(flecsi::exec::accelerator s,
  single<double>::accessor<ro> h_a,
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    field<double>::accessor<ro, na>,
    field<double>::accessor<wo, na>>> v_a,
  std::vector<
    std::tuple<typename vector_field_double<Dim>::template accessor<ro, na>,
      typename vector_field_double<Dim>::template accessor<ro, na>,
      typename vector_field_double<Dim>::template accessor<wo, na>>>
    vec_v_a) noexcept {

  // Scalar
  for(auto & [n_a, dt_a, new_a] : v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, n_a.span().size())) {
      new_a(i) = n_a(i) + *h_a * dt_a(i);
    }; // forall
  }

  // Vec<D>
  for(auto & [n_a, dt_a, new_a] : vec_v_a) {
    auto gpu_n_a = vector_field_double<Dim>::to_accelerator(n_a);
    auto gpu_dt_a = vector_field_double<Dim>::to_accelerator(dt_a);
    auto gpu_new_a = vector_field_double<Dim>::to_accelerator(new_a);
    for(int axis = 0; axis < n_a.size(); axis++) {
      s.executor().forall(
        i, flecsi::util::iota_view({}, gpu_n_a[axis].span().size())) {
        gpu_new_a[axis](i) = gpu_n_a[axis](i) + *h_a * gpu_dt_a[axis](i);
      }; // forall
    }
  }
}

struct update_u_stage_scalar_task : annotation_task_region {
  inline static const std::string name{"update_u_stage_scalar"};
};

template<std::size_t Dim>
void
update_u_stage_scalar(flecsi::exec::accelerator s,
  single<double>::accessor<ro> h_a,
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    field<double>::accessor<ro, na>,
    field<double>::accessor<wo, na>>> v_a) noexcept {

  // Scalar
  for(auto & [n_a, dt_a, new_a] : v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, n_a.span().size())) {
      new_a(i) = n_a(i) + *h_a * dt_a(i);
    }; // forall
  }
}

struct add_k1_k2_task : annotation_task_region {
  inline static const std::string name{"add_k1_k2"};
};

template<std::size_t Dim>
void
add_k1_k2(flecsi::exec::accelerator s,
  std::vector<std::tuple<field<double>::accessor<rw, na>,
    field<double>::accessor<ro, na>>> v_a,
  std::vector<
    std::tuple<typename vector_field_double<Dim>::template accessor<rw, na>,
      typename vector_field_double<Dim>::template accessor<ro, na>>>
    vec_v_a) noexcept {

  // Scalar
  for(auto & [rk_1_a, rk_2_a] : v_a) {
    s.executor().forall(i, flecsi::util::iota_view({}, rk_1_a.span().size())) {
      rk_1_a(i) = (rk_1_a(i) + rk_2_a(i)) * 0.5;
    }; // forall
  }

  // Vec<D>
  for(auto & [rk_1_a, rk_2_a] : vec_v_a) {
    auto gpu_rk_1_a = vector_field_double<Dim>::to_accelerator(rk_1_a);
    auto gpu_rk_2_a = vector_field_double<Dim>::to_accelerator(rk_2_a);
    for(int axis = 0; axis < rk_1_a.size(); axis++) {
      s.executor().forall(
        i, flecsi::util::iota_view({}, gpu_rk_1_a[axis].span().size())) {
        gpu_rk_1_a[axis](i) = (gpu_rk_1_a[axis](i) + gpu_rk_2_a[axis](i)) * 0.5;
      }; // forall
    }
  }
}

} // namespace hard::tasks::hydro

#endif // HARD_MODULE_HYDRO_RHS_HH
