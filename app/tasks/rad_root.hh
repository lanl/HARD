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
  single<double>::accessor<ro> kappa_a,
  single<double>::accessor<ro> particle_mass_a,
  single<double>::accessor<ro> dt_a,
  field<double>::accessor<rw, na> temperature_a,
  field<double>::accessor<rw, na> dt_total_energy_density_implicit_a,
  field<double>::accessor<rw, na> dt_radiation_energy_density_implicit_a,
  eos::eos_wrapper const & eos) {

  using hard::tasks::util::get_mdiota_policy;

  auto & dt_weighted = *dt_a;
  const double one_over_aii_dt{1.0 / dt_weighted};

  auto r = m.template mdcolex<is::cells>(r_a);
  auto u = m.template mdcolex<is::cells>(u_a);
  auto temperature = m.template mdcolex<is::cells>(temperature_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  auto dt_total_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_total_energy_density_implicit_a);
  auto dt_radiation_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_a);

  auto const kappa = *kappa_a;
  auto const particle_mass = *particle_mass_a;

  const double dt_constant =
    hard::constants::cgs::speed_of_light * kappa * dt_weighted;

  const double dt_constant_rc = hard::constants::cgs::speed_of_light * kappa *
                                dt_weighted *
                                hard::constants::cgs::radiation_constant;

  const double one_plus_dt_constant =
    1.0 + hard::constants::cgs::speed_of_light * kappa * dt_weighted;

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
      // auto const gamma = *gamma_a;

      const double ke = 0.5 * r(i) * u(i).norm_squared(); // kinetic energy
      const double en = rE(i) - ke; // internal energy
      
      //getting temperature from EOS, this should be inital guess
      temperature(i) = eos.tRhoSie(r(i), en);

      const double En = radiation_energy_density(i); // radiation energy
      assert(En > 0);
      //Here, we need to get updated Temperature via root-finding
      //We have the form like:
      //F(T) = e(rho,T) - e^n + a/(1+a) * (ar*T^4 - En) 
      // where a = dt*kappa*c
      // For Newton-Raphson method, we need find T such that F(T) == 0
      const double up_Tn = 0.0;
      double TempFour = pow(up_Tn, 4.0);

      const double up_En =
        (En + dt_constant_rc * TempFour) * one_plus_dt_constant;
      // Simply get updated en from T^{n+1} via EOS
      //const double up_en = en + dt_constant_rc * TempFour - dt_constant * up_En;
      const double up_en = eos.eRhoT(r(i), up_Tn);

      dt_total_energy_density_implicit(i) += (up_en - en) * one_over_aii_dt;
      dt_radiation_energy_density_implicit(i) += (up_En - En) * one_over_aii_dt;

    }; // for
  }
  else if constexpr(D == 2) {

    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "upd_energy_density") {

      auto [j, i] = ji;

      const double ke =
        0.5 * r(i, j) * u(i, j).norm_squared(); // kinetic energy
      const double en = rE(i, j) - ke; // internal energy

      temperature(i, j) = eos.tRhoSie(r(i, j), en);
      double TempFour = pow(temperature(i, j), 4.0);

      const double En = radiation_energy_density(i, j); // radiation energy
      assert(En > 0);
      const double up_En =
        (En + dt_constant_rc * TempFour) * one_plus_dt_constant;
      const double up_en = en + dt_constant_rc * TempFour - dt_constant * up_En;

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
      auto [k, j, i] = kji;

      const double ke =
        0.5 * r(i, j, k) * u(i, j, k).norm_squared(); // kinetic energy
      const double en = rE(i, j, k) - ke; // internal energy

      temperature(i, j, k) = eos.tRhoSie(r(i, j, k), en);
      double TempFour = pow(temperature(i, j, k), 4.0);

      const double En = radiation_energy_density(i, j, k); // radiation energy
      assert(En > 0);
      const double up_En =
        (En + dt_constant_rc * TempFour) * one_plus_dt_constant;
      const double up_en = en + dt_constant_rc * TempFour - dt_constant * up_En;

      dt_total_energy_density_implicit(i, j, k) +=
        (up_en - en) * one_over_aii_dt;
      dt_radiation_energy_density_implicit(i, j, k) +=
        (en - up_en) * one_over_aii_dt;
    };
  }
}
} // namespace hard::task::rad_root
#endif
