#ifndef HARD_MODULE_HYDRO_TIME_DERIVATIVE_HH
#define HARD_MODULE_HYDRO_TIME_DERIVATIVE_HH

#include "../modules/hydro/numerical_algorithms/time_stepper.hh"
#include "types.hh"
#include "utils.hh"
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
  typename RK<Dim>::accessor<wo, na> rk_a) noexcept {

  auto [dt_mass, dt_etot, dt_mom] = rk_a;

  s.executor().forall(i, flecsi::util::iota_view({}, dt_mass.span().size())) {
    dt_mass(i) = 0.0;
    dt_mom(i) = vec<Dim>(0.0);
    dt_etot(i) = 0.0;
  }; // forall
}

//
// Store the evolved variables U^n at t=t^n into temporary space before running
// RK substeps.
//
template<std::size_t Dim>
void
store_current_state(flecsi::exec::accelerator s,
  // Copied from
  typename RK<Dim>::accessor<ro, na> rk_1_a,
  typename RK<Dim>::accessor<wo, na> rk_2_a) noexcept {

  auto [mass_1, etot_1, mom_1] = rk_1_a;
  auto [mass_2, etot_2, mom_2] = rk_2_a;

  s.executor().forall(i, flecsi::util::iota_view({}, mass_1.span().size())) {
    mass_2(i) = mass_1(i);
    mom_2(i) = mom_1(i);
    etot_2(i) = etot_1(i);
  }; // forall
}

//
// Used for the RK substeps
//
template<std::size_t Dim>
void
update_u(flecsi::exec::accelerator s,
  single<double>::accessor<ro> dt_a,
  typename mesh<Dim>::template accessor<ro> m,
  // U^n we want to update
  typename RK<Dim>::accessor<rw, na> rk_n_a,
  // Time derivatives for the state U^1
  typename RK<Dim>::accessor<ro, na> rk_dt_a) noexcept {

  auto [mass_density, total_energy_density, momentum_density] =
    RK<Dim>::mdcolex(m, rk_n_a);

  auto [dt_mass_density, dt_total_energy_density, dt_momentum_density] =
    RK<Dim>::mdcolex(m, rk_dt_a);

  using hard::tasks::util::get_mdiota_policy;

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      auto h = *dt_a;
      mass_density(i) += h * dt_mass_density(i);
      momentum_density(i) += h * dt_momentum_density(i);
      total_energy_density(i) += h * dt_total_energy_density(i);

    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto h = *dt_a;
      auto [j, i] = ji;
      // Weights
      mass_density(i, j) += h * dt_mass_density(i, j);
      momentum_density(i, j) += h * dt_momentum_density(i, j);
      total_energy_density(i, j) += h * dt_total_energy_density(i, j);
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto h = *dt_a;
      auto [k, j, i] = kji;

      mass_density(i, j, k) += h * dt_mass_density(i, j, k);
      momentum_density(i, j, k) += h * dt_momentum_density(i, j, k);
      total_energy_density(i, j, k) += h * dt_total_energy_density(i, j, k);
    }; // forall
  }
}

//
// Used for the RK stage updates
//
template<std::size_t Dim, time_stepper::rk_stage stage>
void
update_u_stage(flecsi::exec::cpu s,
  single<double>::accessor<ro> dt_a,
  typename mesh<Dim>::template accessor<ro> m,
  // U^n we want to update
  typename RK<Dim>::accessor<rw, na> rk_n_a,
  // Time derivatives for the state U^1
  typename RK<Dim>::accessor<ro, na> rk_dt1_a,
  // U^n+1 updated after stage
  field<double>::accessor<rw, na> mass_density_b,
  typename field<vec<Dim>>::template accessor<rw, na> momentum_density_b,
  field<double>::accessor<rw, na> total_energy_density_b) noexcept {

  auto [mass_density, total_energy_density, momentum_density] =
    RK<Dim>::mdcolex(m, rk_n_a);

  auto mass_density_new = m.template mdcolex<is::cells>(mass_density_b);
  auto momentum_density_new = m.template mdcolex<is::cells>(momentum_density_b);
  auto total_energy_density_new =
    m.template mdcolex<is::cells>(total_energy_density_b);

  auto [dt_mass_density, dt_total_energy_density, dt_momentum_density] =
    RK<Dim>::mdcolex(m, rk_dt1_a);

  auto h = *dt_a;
  using hard::tasks::util::get_mdiota_policy;

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      mass_density_new(i) = mass_density(i) + h * dt_mass_density(i);
      momentum_density_new(i) =
        momentum_density(i) + h * dt_momentum_density(i);
      total_energy_density_new(i) =
        total_energy_density(i) + h * dt_total_energy_density(i);
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      mass_density_new(i, j) = mass_density(i, j) + h * dt_mass_density(i, j);
      momentum_density_new(i, j) =
        momentum_density(i, j) + h * dt_momentum_density(i, j);
      total_energy_density_new(i, j) =
        total_energy_density(i, j) + h * dt_total_energy_density(i, j);
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      mass_density_new(i, j, k) =
        mass_density(i, j, k) + h * dt_mass_density(i, j, k);
      momentum_density_new(i, j, k) =
        momentum_density(i, j, k) + h * dt_momentum_density(i, j, k);
      total_energy_density_new(i, j, k) =
        total_energy_density(i, j, k) + h * dt_total_energy_density(i, j, k);
    }; // forall
  }
}

template<std::size_t Dim>
void
add_k1_k2(flecsi::exec::accelerator s,
  typename mesh<Dim>::template accessor<ro> m,
  // K1
  typename RK<Dim>::accessor<rw, na> rk_dt1_a,
  // K2
  typename RK<Dim>::accessor<ro, na> rk_dt2_a) noexcept {
  // K1
  auto [dt_r, dt_te, dt_ru] = RK<Dim>::mdcolex(m, rk_dt1_a);

  // K2
  auto [dt_r2, dt_te2, dt_ru2] = RK<Dim>::mdcolex(m, rk_dt2_a);

  using hard::tasks::util::get_mdiota_policy;

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_r(i) = (dt_r(i) + dt_r2(i)) * 0.5;
      dt_ru(i) = (dt_ru(i) + dt_ru2(i)) * 0.5;
      dt_te(i) = (dt_te(i) + dt_te2(i)) * 0.5;
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(dt_r,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      // Weights
      dt_r(i, j) = (dt_r(i, j) + dt_r2(i, j)) * 0.5;
      dt_ru(i, j) = (dt_ru(i, j) + dt_ru2(i, j)) * 0.5;
      dt_te(i, j) = (dt_te(i, j) + dt_te2(i, j)) * 0.5;
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(dt_r,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      dt_r(i, j, k) = (dt_r(i, j, k) + dt_r2(i, j, k)) * 0.5;
      dt_ru(i, j, k) = (dt_ru(i, j, k) + dt_ru2(i, j, k)) * 0.5;
      dt_te(i, j, k) = (dt_te(i, j, k) + dt_te2(i, j, k)) * 0.5;
    }; // forall
  }
}

} // namespace hard::tasks::hydro

#endif // HARD_MODULE_HYDRO_TIME_DERIVATIVE_HH
