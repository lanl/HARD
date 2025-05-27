
#pragma once

#include "../../types.hh"
#include "../utils.hh"
#include "spec/eos.hh"
#include <cstddef>

namespace hard::tasks::hydro {

template<std::size_t Dim>
void
conservative_to_primitive(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<ro, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, na> momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> velocity_a,
  field<double>::accessor<wo, na> pressure_a,
  field<double>::accessor<wo, na> specific_internal_energy_a,
  field<double>::accessor<wo, na> soundspeed_a,
  eos::eos_wrapper const & eos) {
  using hard::tasks::util::get_mdiota_policy;

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
  auto specific_internal_energy =
    m.template mdcolex<is::cells>(specific_internal_energy_a);
  auto soundspeed = m.template mdcolex<is::cells>(soundspeed_a);

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      velocity(i) = momentum_density(i) / mass_density(i);
      specific_internal_energy(i) =
        (total_energy_density(i) -
          0.5 * mass_density(i) * (velocity(i).norm_squared())) /
        mass_density(i);
      pressure(i) = eos.pRhoSie(mass_density(i), specific_internal_energy(i));
      soundspeed(i) = eos.cRhoSie(mass_density(i), specific_internal_energy(i));
    }; // for
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_yx = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_yx) {
      auto [j, i] = ji;
      velocity(i, j) = momentum_density(i, j) / mass_density(i, j);
      specific_internal_energy(i, j) =
        total_energy_density(i, j) / mass_density(i, j) -
        0.5 * velocity(i, j).norm_squared();
      pressure(i, j) =
        eos.pRhoSie(mass_density(i, j), specific_internal_energy(i, j));
      soundspeed(i, j) =
        eos.cRhoSie(mass_density(i, j), specific_internal_energy(i, j));
    }; // for
  }
  else /* Dim == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_zyx) {
      auto [k, j, i] = kji;
      velocity(i, j, k) = momentum_density(i, j, k) / mass_density(i, j, k);
      specific_internal_energy(i, j, k) =
        total_energy_density(i, j, k) / mass_density(i, j, k) -
        0.5 * velocity(i, j, k).norm_squared();
      pressure(i, j, k) =
        eos.pRhoSie(mass_density(i, j, k), specific_internal_energy(i, j, k));
      soundspeed(i, j, k) =
        eos.cRhoSie(mass_density(i, j, k), specific_internal_energy(i, j, k));

    }; // forall
  } // if
}

} // namespace hard::tasks::hydro
