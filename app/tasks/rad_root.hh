#ifndef HARD_TASKS_RADROOT_HH
#define HARD_TASKS_RADROOT_HH

#include "../constants.hh"
#include "../numerical_algorithms/root_finder.hh"
#include "utils.hh"
#include <spec/utils.hh>

namespace hard::task::rad_root {
// This task computes the updated radiation and gas energies.
// The update is performed by finding a root to a quartic polynomial.
template<std::size_t D>
void
update_energy_density(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  typename field<vec<D>>::template accessor<ro, na> u_a,
  field<double>::accessor<ro, na> rE_a,
  field<double>::accessor<ro, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a,
  single<double>::accessor<ro> kappa_a,
  single<double>::accessor<ro> particle_mass_a,
  single<double>::accessor<ro> dt_a,
  //
  field<double>::accessor<rw, na> dt_total_energy_density_implicit_a,
  field<double>::accessor<rw, na> dt_radiation_energy_density_implicit_a) {

  using hard::tasks::util::get_mdiota_policy;

  auto r = m.template mdcolex<is::cells>(r_a);
  auto u = m.template mdcolex<is::cells>(u_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  auto dt_total_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_total_energy_density_implicit_a);
  auto dt_radiation_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_a);

  if constexpr(D == 1) {
    /*------------------------------------------------------------------------*
      Compute the updated gas energy and raditaion energy.
     *------------------------------------------------------------------------*/

    forall(
      i, (m.template cells<ax::x, dm::quantities>()), "upd_energy_density") {
      auto & dt_weighted = *dt_a;
      const double one_over_aii_dt{1.0 / dt_weighted};
      // a1 and a2 are constants defined in the paper Moens et al (2022)

      // TODO: verify these variables are not computed by each thread, and
      // stored in each thread register
      auto const gamma = *gamma_a;
      auto const kappa = *kappa_a;
      auto const particle_mass = *particle_mass_a;
      const double gamma_minus_1 = gamma - 1.0;
      const double a1_times_rho3 =
        4.0 * kappa * hard::constants::cgs::stefan_boltzmann_constant *
        spec::utils::sqr(
          spec::utils::sqr(gamma_minus_1 * particle_mass /
                           hard::constants::cgs::boltzmann_constant)) *
        dt_weighted;
      const double a2_over_rho =
        hard::constants::cgs::speed_of_light * kappa * dt_weighted;

      const double a1 = a1_times_rho3 / (r(i) * spec::utils::sqr((r(i))));
      const double a2 = a2_over_rho * r(i);

      const double ke = 0.5 * r(i) * u(i).norm_squared(); // kinetic energy
      const double en = rE(i) - ke; // internal energy
      const double En = radiation_energy_density(i); // radiation energy

      const double up_en = numerical_algorithms::root_finder::halleys_get_root(
        a1, 1.0 + a2, -(en + (en + En) * a2), en + En);

      dt_total_energy_density_implicit(i) += (up_en - en) * one_over_aii_dt;
      dt_radiation_energy_density_implicit(i) += (en - up_en) * one_over_aii_dt;
    }; // for
  }
  else if constexpr(D == 2) {

    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "upd_energy_density") {
      auto & dt_weighted = *dt_a;
      const double one_over_aii_dt{1.0 / dt_weighted};
      // TODO: verify these variables are not computed by each thread, and
      // stored in each thread register
      auto const gamma = *gamma_a;
      auto const kappa = *kappa_a;
      auto const particle_mass = *particle_mass_a;
      const double gamma_minus_1 = gamma - 1.0;
      const double a1_times_rho3 =
        4.0 * kappa * hard::constants::cgs::stefan_boltzmann_constant *
        spec::utils::sqr(
          spec::utils::sqr(gamma_minus_1 * particle_mass /
                           hard::constants::cgs::boltzmann_constant)) *
        dt_weighted;
      const double a2_over_rho =
        hard::constants::cgs::speed_of_light * kappa * dt_weighted;

      auto [j, i] = ji;
      const double a1 = a1_times_rho3 / (r(i, j) * spec::utils::sqr((r(i, j))));
      const double a2 = a2_over_rho * r(i, j);

      const double ke = 0.5 * r(i, j) * u(i, j).norm_squared();
      const double en = rE(i, j) - ke; // internal energy
      const double En = radiation_energy_density(i, j); // radiation energy

      const double up_en = numerical_algorithms::root_finder::halleys_get_root(
        a1, 1.0 + a2, -(en + (en + En) * a2), en + En);

      dt_total_energy_density_implicit(i, j) += (up_en - en) * one_over_aii_dt;
      dt_radiation_energy_density_implicit(i, j) +=
        (en - up_en) * one_over_aii_dt;
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "upd_energy_density") {
      auto & dt_weighted = *dt_a;
      const double one_over_aii_dt{1.0 / dt_weighted};
      // TODO: verify these variables are not computed by each thread, and
      // stored in each thread register
      auto const gamma = *gamma_a;
      auto const kappa = *kappa_a;
      auto const particle_mass = *particle_mass_a;
      const double gamma_minus_1 = gamma - 1.0;
      const double a1_times_rho3 =
        4.0 * kappa * hard::constants::cgs::stefan_boltzmann_constant *
        spec::utils::sqr(
          spec::utils::sqr(gamma_minus_1 * particle_mass /
                           hard::constants::cgs::boltzmann_constant)) *
        dt_weighted;
      const double a2_over_rho =
        hard::constants::cgs::speed_of_light * kappa * dt_weighted;

      auto [k, j, i] = kji;
      const double a1 =
        a1_times_rho3 / (r(i, j, k) * spec::utils::sqr((r(i, j, k))));
      const double a2 = a2_over_rho * r(i, j, k);

      const double ke = 0.5 * r(i, j, k) * u(i, j, k).norm_squared();
      const double en = rE(i, j, k) - ke; // internal energy
      const double En = radiation_energy_density(i, j, k); // radiation energy

      const double up_en = numerical_algorithms::root_finder::halleys_get_root(
        a1, 1.0 + a2, -(en + (en + En) * a2), en + En);

      dt_total_energy_density_implicit(i, j, k) +=
        (up_en - en) * one_over_aii_dt;
      dt_radiation_energy_density_implicit(i, j, k) +=
        (en - up_en) * one_over_aii_dt;
    };
  }
}
} // namespace hard::task::rad_root
#endif
