
#ifndef HARD_TIME_DERIVATIVE_HH
#define HARD_TIME_DERIVATIVE_HH

#include "../numerical_algorithms/time_stepper.hh"
#include "../types.hh"
#include "utils.hh"
#include <cstddef>

namespace hard::tasks {

template<std::size_t Dim>
void
set_dudt_to_zero(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> dt_mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> dt_momentum_density_a,
  field<double>::accessor<wo, na> dt_total_energy_density_a,
  field<double>::accessor<wo, na> dt_radiation_energy_density_a,
  field<double>::accessor<wo, na> dt_total_energy_density_implicit_a,
  field<double>::accessor<wo, na> dt_radiation_energy_density_implicit_a) {

  auto dt_mass_density = m.template mdcolex<is::cells>(dt_mass_density_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);
  auto dt_total_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_total_energy_density_implicit_a);
  auto dt_radiation_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_a);

  using hard::tasks::util::get_mdiota_policy;
  if constexpr(Dim == 1) {
    forall(
      i, (m.template cells<ax::x, dm::quantities>()), "set_dudt_to_zero_1d") {
      dt_mass_density(i) = 0.0;
      dt_momentum_density(i).x = 0.0;
      dt_total_energy_density(i) = 0.0;
      dt_radiation_energy_density(i) = 0.0;
      dt_total_energy_density_implicit(i) = 0.0;
      dt_radiation_energy_density_implicit(i) = 0.0;
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "set_dudt_to_zero_2d") {
      auto [j, i] = ji;
      dt_mass_density(i, j) = 0.0;
      dt_momentum_density(i, j).x = 0.0;
      dt_momentum_density(i, j).y = 0.0;
      dt_total_energy_density(i, j) = 0.0;
#ifndef DISABLE_RADIATION
      dt_radiation_energy_density(i, j) = 0.0;
      dt_total_energy_density_implicit(i, j) = 0.0;
      dt_radiation_energy_density_implicit(i, j) = 0.0;
#endif

    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(dt_mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "set_dudt_to_zero_3d") {
      auto [k, j, i] = kji;
      dt_mass_density(i, j, k) = 0.0;
      dt_momentum_density(i, j, k).x = 0.0;
      dt_momentum_density(i, j, k).y = 0.0;
      dt_momentum_density(i, j, k).z = 0.0;
      dt_total_energy_density(i, j, k) = 0.0;
#ifndef DISABLE_RADIATION
      dt_radiation_energy_density(i, j, k) = 0.0;
      dt_total_energy_density_implicit(i, j, k) = 0.0;
      dt_radiation_energy_density_implicit(i, j, k) = 0.0;
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
store_current_state(typename mesh<Dim>::template accessor<ro> m,
  // Copied from
  field<double>::accessor<ro, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, na> momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    radiation_energy_density_a
#endif
  // Copied into
  ,
  field<double>::accessor<wo, na> mass_density_n_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_n_a,
  field<double>::accessor<wo, na> total_energy_density_n_a,
  field<double>::accessor<wo, na>
#ifndef DISABLE_RADIATION
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
#ifndef DISABLE_RADIATION
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto radiation_energy_density_n =
    m.template mdcolex<is::cells>(radiation_energy_density_n_a);
#endif

  using hard::tasks::util::get_mdiota_policy;
  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "store_state_1d") {
      mass_density_n(i) = mass_density(i);
      momentum_density_n(i) = momentum_density(i);
      total_energy_density_n(i) = total_energy_density(i);
#ifndef DISABLE_RADIATION
      radiation_energy_density_n(i) = radiation_energy_density(i);
#endif
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "store_state_2d") {
      auto [j, i] = ji;
      mass_density_n(i, j) = mass_density(i, j);
      momentum_density_n(i, j) = momentum_density(i, j);
      total_energy_density_n(i, j) = total_energy_density(i, j);
#ifndef DISABLE_RADIATION
      radiation_energy_density_n(i, j) = radiation_energy_density(i, j);
#endif

    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "store_state_3d") {
      auto [k, j, i] = kji;
      mass_density_n(i, j, k) = mass_density(i, j, k);
      momentum_density_n(i, j, k) = momentum_density(i, j, k);
      total_energy_density_n(i, j, k) = total_energy_density(i, j, k);
#ifndef DISABLE_RADIATION
      radiation_energy_density_n(i, j, k) = radiation_energy_density(i, j, k);
#endif

    }; // forall
  }
}

//
// Prepare evolved variables to U^1 or U^2 after implicit steps (radiative
// heating / diffusion) and storing (dU_dt)_implicit have been done.
//
template<std::size_t Dim>
void
compute_u_after_implicit_solve(typename mesh<Dim>::template accessor<ro> m,
  single<double>::accessor<ro> dtw_a,
  field<double>::accessor<rw, na> total_energy_density_a,
  field<double>::accessor<rw, na>
#ifndef DISABLE_RADIATION
    radiation_energy_density_a
#endif
  ,
  field<double>::accessor<ro, na> dt_total_energy_density_implicit_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    dt_radiation_energy_density_implicit_a
#endif
) {

  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
#ifndef DISABLE_RADIATION

  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
#endif
  auto dt_total_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_total_energy_density_implicit_a);
#ifndef DISABLE_RADIATION
  auto dt_radiation_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_a);
#endif

  using hard::tasks::util::get_mdiota_policy;
  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "compute_u_1d") {
      auto & dt_weighted = *dtw_a;
      total_energy_density(i) +=
        dt_weighted * dt_total_energy_density_implicit(i);
#ifndef DISABLE_RADIATION
      radiation_energy_density(i) +=
        dt_weighted * dt_radiation_energy_density_implicit(i);
#endif
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(total_energy_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "compute_u_2d") {
      auto [j, i] = ji;
      auto & dt_weighted = *dtw_a;
      total_energy_density(i, j) +=
        dt_weighted * dt_total_energy_density_implicit(i, j);
#ifndef DISABLE_RADIATION
      radiation_energy_density(i, j) +=
        dt_weighted * dt_radiation_energy_density_implicit(i, j);
#endif
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(total_energy_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "compute_u_3d") {
      auto [k, j, i] = kji;
      auto & dt_weighted = *dtw_a;
      total_energy_density(i, j, k) +=
        dt_weighted * dt_total_energy_density_implicit(i, j, k);
#ifndef DISABLE_RADIATION
      radiation_energy_density(i, j, k) +=
        dt_weighted * dt_radiation_energy_density_implicit(i, j, k);
#endif
    }; // forall
  }
}

//
// Used for the RK substep 2 and the final update
//
//  - If (Stage == Second), prepare an intermediate state
//      U^n + E(U^1) + (1 - 2g) * I(U^1)
//
//  - If (Stage == Update), perform the time update
//      U^{n+1} = U^n + 0.5 * (E(U^1) + E(U^2) + I(U^1) + I(U^2))
//
template<std::size_t Dim, time_stepper::rk_stage Stage>
void
update_u(single<double>::accessor<ro> dt_a,
  typename mesh<Dim>::template accessor<ro> m,
  // U^n
  field<double>::accessor<ro, na> mass_density_n_a,
  typename field<vec<Dim>>::template accessor<ro, na> momentum_density_n_a,
  field<double>::accessor<ro, na> total_energy_density_n_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    radiation_energy_density_n_a
#endif
  ,
  // explicit & implicit time derivatives for the state U^1
  field<double>::accessor<ro, na> dt_mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, na> dt_momentum_density_a,
  field<double>::accessor<ro, na> dt_total_energy_density_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    dt_radiation_energy_density_a
#endif
  ,
  field<double>::accessor<ro, na> dt_total_energy_density_implicit_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    dt_radiation_energy_density_implicit_a
#endif
  ,
  // explicit & implicit time derivatives for the state U^2
  field<double>::accessor<ro, na> dt_mass_density_2_a,
  typename field<vec<Dim>>::template accessor<ro, na> dt_momentum_density_2_a,
  field<double>::accessor<ro, na> dt_total_energy_density_2_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    dt_radiation_energy_density_2_a
#endif
  ,
  field<double>::accessor<ro, na> dt_total_energy_density_implicit_2_a,
  field<double>::accessor<ro, na>
#ifndef DISABLE_RADIATION
    dt_radiation_energy_density_implicit_2_a
#endif
  ,
  // write U^n+1 into following arguments
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na>
#ifndef DISABLE_RADIATION
    radiation_energy_density_a
#endif
) {

  auto mass_density_n = m.template mdcolex<is::cells>(mass_density_n_a);
  auto momentum_density_n = m.template mdcolex<is::cells>(momentum_density_n_a);
  auto total_energy_density_n =
    m.template mdcolex<is::cells>(total_energy_density_n_a);
#ifndef DISABLE_RADIATION
  auto radiation_energy_density_n =
    m.template mdcolex<is::cells>(radiation_energy_density_n_a);
#endif

  auto dt_mass_density = m.template mdcolex<is::cells>(dt_mass_density_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
#ifndef DISABLE_RADIATION
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);
#endif
  auto dt_total_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_total_energy_density_implicit_a);
#ifndef DISABLE_RADIATION
  auto dt_radiation_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_a);
#endif

  auto dt_mass_density_2 = m.template mdcolex<is::cells>(dt_mass_density_2_a);
  auto dt_momentum_density_2 =
    m.template mdcolex<is::cells>(dt_momentum_density_2_a);
  auto dt_total_energy_density_2 =
    m.template mdcolex<is::cells>(dt_total_energy_density_2_a);
#ifndef DISABLE_RADIATION
  auto dt_radiation_energy_density_2 =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_2_a);
#endif
  auto dt_total_energy_density_implicit_2 =
    m.template mdcolex<is::cells>(dt_total_energy_density_implicit_2_a);
#ifndef DISABLE_RADIATION
  auto dt_radiation_energy_density_implicit_2 =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_2_a);
#endif

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
#ifndef DISABLE_RADIATION
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
#endif
  using hard::tasks::util::get_mdiota_policy;
  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "update_u_1d") {
      // Weights
      auto dt = *dt_a;
      const double one_minus_2gamma_dt =
        time_stepper::one_minus_2_time_stepper_gamma * dt;

      if constexpr(Stage == time_stepper::rk_stage::Second) {

        mass_density(i) = mass_density_n(i) + (dt * dt_mass_density(i));

        momentum_density(i) =
          momentum_density_n(i) + (dt * dt_momentum_density(i));

        total_energy_density(i) =
          total_energy_density_n(i) + (dt * dt_total_energy_density(i)) +
          (one_minus_2gamma_dt * dt_total_energy_density_implicit(i));

#ifndef DISABLE_RADIATION
        radiation_energy_density(i) =
          radiation_energy_density_n(i) +
          (dt * dt_radiation_energy_density(i)) +
          (one_minus_2gamma_dt * dt_radiation_energy_density_implicit(i));
#endif
      }
      else if constexpr(Stage == time_stepper::rk_stage::Update) {

        mass_density(i) =
          mass_density_n(i) +
          0.5 * dt * (dt_mass_density(i) + dt_mass_density_2(i));

        momentum_density(i) =
          momentum_density_n(i) +
          0.5 * dt * (dt_momentum_density(i) + dt_momentum_density_2(i));

        total_energy_density(i) =
          total_energy_density_n(i) +
          0.5 * dt *
            (dt_total_energy_density(i) + dt_total_energy_density_2(i) +
              dt_total_energy_density_implicit(i) +
              dt_total_energy_density_implicit_2(i));

#ifndef DISABLE_RADIATION
        radiation_energy_density(i) =
          radiation_energy_density_n(i) +
          0.5 * dt *
            (dt_radiation_energy_density(i) + dt_radiation_energy_density_2(i) +
              dt_radiation_energy_density_implicit(i) +
              dt_radiation_energy_density_implicit_2(i));
#endif
      }
      else {
        flog_fatal("Incorrect call of update_u()");
      }
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_qq = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "update_u_2d") {
      auto [j, i] = ji;
      // Weights
      auto dt = *dt_a;
      const double one_minus_2gamma_dt =
        time_stepper::one_minus_2_time_stepper_gamma * dt;

      if constexpr(Stage == time_stepper::rk_stage::Second) {

        mass_density(i, j) =
          mass_density_n(i, j) + (dt * dt_mass_density(i, j));

        momentum_density(i, j) =
          momentum_density_n(i, j) + (dt * dt_momentum_density(i, j));

        total_energy_density(i, j) =
          total_energy_density_n(i, j) + (dt * dt_total_energy_density(i, j)) +
          (one_minus_2gamma_dt * dt_total_energy_density_implicit(i, j));

#ifndef DISABLE_RADIATION
        radiation_energy_density(i, j) =
          radiation_energy_density_n(i, j) +
          (dt * dt_radiation_energy_density(i, j)) +
          (one_minus_2gamma_dt * dt_radiation_energy_density_implicit(i, j));
#endif
      }
      else if constexpr(Stage == time_stepper::rk_stage::Update) {

        mass_density(i, j) =
          mass_density_n(i, j) +
          0.5 * dt * (dt_mass_density(i, j) + dt_mass_density_2(i, j));

        momentum_density(i, j) =
          momentum_density_n(i, j) +
          0.5 * dt * (dt_momentum_density(i, j) + dt_momentum_density_2(i, j));

        total_energy_density(i, j) =
          total_energy_density_n(i, j) +
          0.5 * dt *
            (dt_total_energy_density(i, j) + dt_total_energy_density_2(i, j) +
              dt_total_energy_density_implicit(i, j) +
              dt_total_energy_density_implicit_2(i, j));

#ifndef DISABLE_RADIATION
        radiation_energy_density(i, j) =
          radiation_energy_density_n(i, j) +
          0.5 * dt *
            (dt_radiation_energy_density(i, j) +
              dt_radiation_energy_density_2(i, j) +
              dt_radiation_energy_density_implicit(i, j) +
              dt_radiation_energy_density_implicit_2(i, j));
#endif
      }
      else {
        flog_fatal("Incorrect call of update_u()");
      }
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "update_u_3d") {
      auto [k, j, i] = kji;
      // Weights
      auto dt = *dt_a;
      const double one_minus_2gamma_dt =
        time_stepper::one_minus_2_time_stepper_gamma * dt;

      if constexpr(Stage == time_stepper::rk_stage::Second) {
        mass_density(i, j, k) =
          mass_density_n(i, j, k) + (dt * dt_mass_density(i, j, k));

        momentum_density(i, j, k) =
          momentum_density_n(i, j, k) + (dt * dt_momentum_density(i, j, k));

        total_energy_density(i, j, k) =
          total_energy_density_n(i, j, k) +
          (dt * dt_total_energy_density(i, j, k)) +
          (one_minus_2gamma_dt * dt_total_energy_density_implicit(i, j, k));

#ifndef DISABLE_RADIATION
        radiation_energy_density(i, j, k) =
          radiation_energy_density_n(i, j, k) +
          (dt * dt_radiation_energy_density(i, j, k)) +
          (one_minus_2gamma_dt * dt_radiation_energy_density_implicit(i, j, k));
#endif
      }
      else if constexpr(Stage == time_stepper::rk_stage::Update) {
        mass_density(i, j, k) =
          mass_density_n(i, j, k) +
          0.5 * dt * (dt_mass_density(i, j, k) + dt_mass_density_2(i, j, k));

        momentum_density(i, j, k) =
          momentum_density_n(i, j, k) +
          0.5 * dt *
            (dt_momentum_density(i, j, k) + dt_momentum_density_2(i, j, k));

        total_energy_density(i, j, k) =
          total_energy_density_n(i, j, k) +
          0.5 * dt *
            (dt_total_energy_density(i, j, k) +
              dt_total_energy_density_2(i, j, k) +
              dt_total_energy_density_implicit(i, j, k) +
              dt_total_energy_density_implicit_2(i, j, k));

#ifndef DISABLE_RADIATION
        radiation_energy_density(i, j, k) =
          radiation_energy_density_n(i, j, k) +
          0.5 * dt *
            (dt_radiation_energy_density(i, j, k) +
              dt_radiation_energy_density_2(i, j, k) +
              dt_radiation_energy_density_implicit(i, j, k) +
              dt_radiation_energy_density_implicit_2(i, j, k));
#endif
      }
      else {
        flog_fatal("Incorrect call of update_u()");
      }
    }; // forall
  }
}

} // namespace hard::tasks

#endif // HARD_TIME_DERIVATIVE_HH
