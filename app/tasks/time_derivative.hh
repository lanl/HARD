
#ifndef HARD_TIME_DERIVATIVE_HH
#define HARD_TIME_DERIVATIVE_HH

#include "../numerical_algorithms/time_stepper.hh"
#include "../types.hh"
#include "utils.hh"
#include <cstddef>

namespace hard::tasks {

template<std::size_t Dim>
void
set_dudt_to_zero(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> dt_mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> dt_momentum_density_a,
  field<double>::accessor<wo, na> dt_total_energy_density_a,
  field<double>::accessor<wo, na> dt_radiation_energy_density_a) {

  auto dt_mass_density = m.template mdcolex<is::cells>(dt_mass_density_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);

  using hard::tasks::util::get_mdiota_policy;
  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_mass_density(i) = 0.0;
      dt_momentum_density(i).x = 0.0;
      dt_total_energy_density(i) = 0.0;
      dt_radiation_energy_density(i) = 0.0;
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      dt_mass_density(i, j) = 0.0;
      dt_momentum_density(i, j).x = 0.0;
      dt_momentum_density(i, j).y = 0.0;
      dt_total_energy_density(i, j) = 0.0;
#ifdef ENABLE_RADIATION
      dt_radiation_energy_density(i, j) = 0.0;
#endif

    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      dt_mass_density(i, j, k) = 0.0;
      dt_momentum_density(i, j, k).x = 0.0;
      dt_momentum_density(i, j, k).y = 0.0;
      dt_momentum_density(i, j, k).z = 0.0;
      dt_total_energy_density(i, j, k) = 0.0;
#ifdef ENABLE_RADIATION
      dt_radiation_energy_density(i, j, k) = 0.0;
#endif
    }; // forall
  }
}

//
// Store the evolved variables U^n at t=t^n into temporary space before running
// RK substeps.
//
template<std::size_t Dim>
void
store_current_state(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  // Copied from
  field<double>::accessor<ro, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, na> momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  field<double>::accessor<ro, na>
#ifdef ENABLE_RADIATION
    radiation_energy_density_a
#endif
  // Copied into
  ,
  field<double>::accessor<wo, na> mass_density_n_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_n_a,
  field<double>::accessor<wo, na> total_energy_density_n_a,
  field<double>::accessor<wo, na>
#ifdef ENABLE_RADIATION
    radiation_energy_density_n_a
#endif
) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto mass_density_n = m.template mdcolex<is::cells>(mass_density_n_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto momentum_density_n = m.template mdcolex<is::cells>(momentum_density_n_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto total_energy_density_n =
    m.template mdcolex<is::cells>(total_energy_density_n_a);
#ifdef ENABLE_RADIATION
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto radiation_energy_density_n =
    m.template mdcolex<is::cells>(radiation_energy_density_n_a);
#endif

  using hard::tasks::util::get_mdiota_policy;
  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      mass_density_n(i) = mass_density(i);
      momentum_density_n(i) = momentum_density(i);
      total_energy_density_n(i) = total_energy_density(i);
#ifdef ENABLE_RADIATION
      radiation_energy_density_n(i) = radiation_energy_density(i);
#endif
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      mass_density_n(i, j) = mass_density(i, j);
      momentum_density_n(i, j) = momentum_density(i, j);
      total_energy_density_n(i, j) = total_energy_density(i, j);
#ifdef ENABLE_RADIATION
      radiation_energy_density_n(i, j) = radiation_energy_density(i, j);
#endif

    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      mass_density_n(i, j, k) = mass_density(i, j, k);
      momentum_density_n(i, j, k) = momentum_density(i, j, k);
      total_energy_density_n(i, j, k) = total_energy_density(i, j, k);
#ifdef ENABLE_RADIATION
      radiation_energy_density_n(i, j, k) = radiation_energy_density(i, j, k);
#endif

    }; // forall
  }
}

//
// Used for the RK substeps
//
template<std::size_t Dim>
void
update_u(flecsi::exec::cpu s,
  single<double>::accessor<ro> dt_a,
  typename mesh<Dim>::template accessor<ro> m,
  // U^n we want to update
  field<double>::accessor<rw, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, na> momentum_density_a,
  field<double>::accessor<rw, na> total_energy_density_a,
  field<double>::accessor<rw, na>
#ifdef ENABLE_RADIATION
    radiation_energy_density_a
#endif
  ,
  // Time derivatives for the state U^1
  field<double>::accessor<ro, na> dt_mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, na> dt_momentum_density_a,
  field<double>::accessor<ro, na> dt_total_energy_density_a,
  field<double>::accessor<ro, na>
#ifdef ENABLE_RADIATION
    dt_radiation_energy_density_a
#endif
) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
#ifdef ENABLE_RADIATION
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
#endif

  auto dt_mass_density = m.template mdcolex<is::cells>(dt_mass_density_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
#ifdef ENABLE_RADIATION
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);
#endif

  auto h = *dt_a;
  using hard::tasks::util::get_mdiota_policy;

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      mass_density(i) += h * dt_mass_density(i);
      momentum_density(i) += h * dt_momentum_density(i);
      total_energy_density(i) += h * dt_total_energy_density(i);

#ifdef ENABLE_RADIATION
      radiation_energy_density(i) += h * dt_radiation_energy_density(i);
#endif
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      // Weights
      mass_density(i, j) += h * dt_mass_density(i, j);
      momentum_density(i, j) += h * dt_momentum_density(i, j);
      total_energy_density(i, j) += h * dt_total_energy_density(i, j);
#ifdef ENABLE_RADIATION
      radiation_energy_density(i, j) += h * dt_radiation_energy_density(i, j);
#endif
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      mass_density(i, j, k) += h * dt_mass_density(i, j, k);
      momentum_density(i, j, k) += h * dt_momentum_density(i, j, k);
      total_energy_density(i, j, k) += h * dt_total_energy_density(i, j, k);
#ifdef ENABLE_RADIATION
      radiation_energy_density(i, j, k) +=
        h * dt_radiation_energy_density(i, j, k);
#endif
    }; // forall
  }
}

template<std::size_t Dim>
void
add_k1_k2(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  // K1
  field<double>::accessor<rw, na> dt_mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, na> dt_momentum_density_a,
  field<double>::accessor<rw, na> dt_total_energy_density_a,
  field<double>::accessor<rw, na>
#ifdef ENABLE_RADIATION
    dt_radiation_energy_density_a
#endif
  ,
  // K2
  field<double>::accessor<ro, na> dt_mass_density_2_a,
  typename field<vec<Dim>>::template accessor<ro, na> dt_momentum_density_2_a,
  field<double>::accessor<ro, na> dt_total_energy_density_2_a,
  field<double>::accessor<ro, na>
#ifdef ENABLE_RADIATION
    dt_radiation_energy_density_2_a
#endif
) {
  // K1
  auto dt_r = m.template mdcolex<is::cells>(dt_mass_density_a);
  auto dt_ru = m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_te = m.template mdcolex<is::cells>(dt_total_energy_density_a);
#ifdef ENABLE_RADIATION
  auto dt_re = m.template mdcolex<is::cells>(dt_radiation_energy_density_a);
#endif

  // K2
  auto dt_r2 = m.template mdcolex<is::cells>(dt_mass_density_2_a);
  auto dt_ru2 = m.template mdcolex<is::cells>(dt_momentum_density_2_a);
  auto dt_te2 = m.template mdcolex<is::cells>(dt_total_energy_density_2_a);
#ifdef ENABLE_RADIATION
  auto dt_re2 = m.template mdcolex<is::cells>(dt_radiation_energy_density_2_a);
#endif

  using hard::tasks::util::get_mdiota_policy;

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      dt_r(i) = (dt_r(i) + dt_r2(i)) * 0.5;
      dt_ru(i) = (dt_ru(i) + dt_ru2(i)) * 0.5;
      dt_te(i) = (dt_te(i) + dt_te2(i)) * 0.5;

#ifdef ENABLE_RADIATION
      dt_re(i) = (dt_re(i) + dt_re2(i)) * 0.5;
#endif
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
      dt_te(i, j) = (dt_te(i, j) + dt_te2(i, j) * 0.5);
#ifdef ENABLE_RADIATION
      dt_re(i, j) += dt_re2(i, j);
#endif
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
#ifdef ENABLE_RADIATION
      dt_re(i, j, k) = (dt_re(i, j, k) + dt_re2(i, j, k)) * 0.5;
#endif
    }; // forall
  }
}

} // namespace hard::tasks

#endif // HARD_TIME_DERIVATIVE_HH
