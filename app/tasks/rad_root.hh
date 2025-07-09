#ifndef HARD_TASKS_RADROOT_HH
#define HARD_TASKS_RADROOT_HH

#include "../constants.hh"
#include "utils.hh"
#include <spec/utils.hh>

namespace hard::task::rad_root {
// This task computes the updated radiation and gas energies.
// The update is performed by finding a root to a quartic polynomial.

//
template<std::size_t D>
void
update_energy_density(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  typename field<vec<D>>::template accessor<ro, na> u_a,
  field<double>::accessor<rw, na> temperature_a,
  field<double>::accessor<rw, na> rE_a,
  field<double>::accessor<rw, na> radiation_energy_density_a,
  single<double>::accessor<ro> kappa_a,
  single<double>::accessor<ro> dt_a,
  eos::eos_wrapper const & eos) noexcept {

  using hard::tasks::util::get_mdiota_policy;

  auto r = m.template mdcolex<is::cells>(r_a);
  auto u = m.template mdcolex<is::cells>(u_a);
  auto temperature = m.template mdcolex<is::cells>(temperature_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto const kappa = *kappa_a;
  auto & dt_weighted = *dt_a;

  // Timestep constant for energy update
  const double dt_c{hard::constants::cgs::speed_of_light * kappa * dt_weighted};

  // Radiation constant (shorter name)
  const double rad_c{hard::constants::cgs::radiation_constant};

  // Define the En update lambda
  auto get_up_En = [eos, dt_c, rad_c](double_t t, double_t En) {
    return (En + dt_c * rad_c * pow(t, 4.0)) / (1 + dt_c);
  };

  if constexpr(D == 1) {
    /*------------------------------------------------------------------------*
      Compute the updated gas energy and raditaion energy.
     *------------------------------------------------------------------------*/

    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      auto & dt_weighted = *dt_a;
      // NOTE: a1 and a2 are constants defined in the paper Moens et al (2022)
      // FIXME: The variables "a1" and "a2" do not exist, remove comment
      // above?

      // TODO: verify these variables are not computed by each thread, and
      // stored in each thread register
      // auto const gamma = *gamma_a;

      const double ke = 0.5 * r(i) * u(i).norm_squared(); // kinetic energy
      const double en = rE(i) - ke; // internal energy
      assert(en >= 0 && "Internal energy is negative");

      // getting temperature from EOS, this should be inital guess
      temperature(i) = eos.tRhoSie(r(i), en);

      const double En = radiation_energy_density(i); // radiation energy
      assert(En >= 0 && "Radiation energy is negative");

      // Find the next temperature with root finding
      const double up_Tn{tasks::util::find_temp(
        eos, r(i), en, temperature(i), kappa, En, dt_weighted)};

      // Update En with the formula and en with conservation of energy
      const double up_En{get_up_En(up_Tn, En)};
      assert(up_En >= 0 && "Updated radiation energy is negative");

      const double up_en{en - (up_En - En)};
      assert(up_en >= 0 && "Updated internal energy is negative");

      radiation_energy_density(i) = up_En;
      rE(i) = up_en + ke;
      temperature(i) = up_Tn;
    }; // for
  }
  else if constexpr(D == 2) {

    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {

      auto [j, i] = ji;

      const double ke =
        0.5 * r(i, j) * u(i, j).norm_squared(); // kinetic energy
      const double en = rE(i, j) - ke; // internal energy
      assert(en >= 0 && "Internal energy is negative");

      temperature(i, j) = eos.tRhoSie(r(i, j), en);

      const double En = radiation_energy_density(i, j); // radiation energy
      assert(En >= 0 && "Radiation energy is negative");

      // Find the next temperature with root finding
      const double up_Tn{tasks::util::find_temp(
        eos, r(i, j), en, temperature(i, j), kappa, En, dt_weighted)};

      const double up_En{get_up_En(up_Tn, En)};
      assert(up_En >= 0 && "Updated radiation energy is negative");

      const double up_en{en - (up_En - En)};
      assert(up_en >= 0 && "Updated internal energy is negative");

      radiation_energy_density(i, j) = up_En;
      rE(i, j) = up_en + ke;
      temperature(i, j) = up_Tn;
    }; // forall
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      const double ke =
        0.5 * r(i, j, k) * u(i, j, k).norm_squared(); // kinetic energy
      const double en = rE(i, j, k) - ke; // internal energy
      assert(en >= 0 && "Internal energy is negative");

      temperature(i, j, k) = eos.tRhoSie(r(i, j, k), en);

      const double En = radiation_energy_density(i, j, k); // radiation energy
      assert(En >= 0 && "Radiation energy is negative");

      // Find the next temperature with root finding
      const double up_Tn{tasks::util::find_temp(
        eos, r(i, j, k), en, temperature(i, j, k), kappa, En, dt_weighted)};

      const double up_En{get_up_En(up_Tn, En)};
      assert(up_En >= 0 && "Updated radiation energy is negative");

      const double up_en{en - (up_En - En)};
      assert(up_en >= 0 && "Updated internal energy is negative");

      radiation_energy_density(i, j, k) = up_En;
      rE(i, j, k) = up_en + ke;
      temperature(i, j, k) = up_Tn;
    };
  }
}
} // namespace hard::task::rad_root
#endif
