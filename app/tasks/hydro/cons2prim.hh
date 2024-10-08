
#pragma once

#include "../../types.hh"
#include "../utils.hh"
#include <cstddef>

namespace hard::tasks::hydro {

template<std::size_t Dim>
void
conservative_to_primitive(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<ro, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, na> momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> velocity_a,
  field<double>::accessor<wo, na> pressure_a,
  single<double>::accessor<ro> gamma_a) {
  using hard::tasks::util::get_mdiota_policy;

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);

  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "conv2prim") {
      auto const gamma = *gamma_a;
      velocity(i) = momentum_density(i) / mass_density(i);
      pressure(i) =
        (gamma - 1.0) * (total_energy_density(i) -
                          0.5 * mass_density(i) * velocity(i).norm_squared());
    }; // for
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_yx = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_yx, "conv2prim") {
      auto const gamma = *gamma_a;
      auto [j, i] = ji;
      velocity(i, j) = momentum_density(i, j) / mass_density(i, j);
      pressure(i, j) = (gamma - 1.0) * (total_energy_density(i, j) -
                                         0.5 * mass_density(i, j) *
                                           velocity(i, j).norm_squared());
    }; // for
  }
  else /* Dim == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_zyx, "conv2prim") {
      auto const gamma = *gamma_a;
      auto [k, j, i] = kji;
      velocity(i, j, k) = momentum_density(i, j, k) / mass_density(i, j, k);
      pressure(i, j, k) = (gamma - 1.0) * (total_energy_density(i, j, k) -
                                            0.5 * mass_density(i, j, k) *
                                              velocity(i, j, k).norm_squared());
    }; // forall
  } // if
}

} // namespace hard::tasks::hydro
