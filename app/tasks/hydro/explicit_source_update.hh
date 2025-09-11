
#pragma once

#include "../../numerical_algorithms/riemann_solvers.hh"
#include "../../types.hh"
#include "../utils.hh"
#include <cstddef>

namespace hard::tasks::hydro {

using hard::tasks::util::get_mdiota_policy;

// Gravity force and work terms. (Explicit source terms) (for RT instability
// test case)
template<std::size_t D>
void
explicitSourceUpdate(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  // Primitive variables
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  // Terms required for computing radiation force and photon tiring
  typename field<vec<D>>::template accessor<ro, na> gravity_force_a,
  // time derivative
  typename field<vec<D>>::template accessor<rw, na> dt_momentum_density_a,
  field<double>::accessor<rw, na> dt_total_energy_density_a) noexcept {

  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto fg = m.template mdcolex<is::cells>(gravity_force_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
  // const double radiation_constant = hard::constants::cgs::radiation_constant;

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {

      // Explicitly updating conservative variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i) += fg(i);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i) += velocity(i).x() * fg(i).x();

      // TODO:
      // Add the source from the temperature
      // NOTE: Isn't this already in the rad_root part?
      // dt_radiation_energy_density(i) +=
      //   constants::cgs::radiation_constant * pow(T_source(i), 4);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      // Explicitly updating conservative variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j) += fg(i, j);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i, j) +=
        velocity(i, j).x() * fg(i, j).x() + velocity(i, j).y() * fg(i, j).y();
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      // Explicitly updating conservative variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j, k) += fg(i, j, k);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i, j, k) +=
        velocity(i, j, k).x() * fg(i, j, k).x() +
        velocity(i, j, k).y() * fg(i, j, k).y() +
        velocity(i, j, k).z() * fg(i, j, k).z();
    }; // forall
  }
} // explicitSourceUpdate

} // namespace hard::tasks::hydro
