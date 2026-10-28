#ifndef HARD_MODULE_HYDRO_EXTERNAL_SOURCE_HH
#define HARD_MODULE_HYDRO_EXTERNAL_SOURCE_HH

#include "../modules/common/tasks/utils.hh"

namespace hard::tasks {

struct external_source_task : annotation_task_region {
  inline static const std::string name{"external_source"};
};

// Gravity force and work terms. (Explicit source terms) (for RT instability
// test case)
template<std::size_t D>
void
external_source(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  // Primitive variables
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  // Terms required for computing radiation force and photon tiring
  field<double>::accessor<rw, na> mass_density_a,
  typename field<vec<D>>::template accessor<ro, na> gravity_accel_a,
  // time derivative
  typename field<vec<D>>::template accessor<rw, na> dt_momentum_density_a,
  field<double>::accessor<rw, na> dt_total_energy_density_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto ga = m.template mdcolex<is::cells>(gravity_accel_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);

  if constexpr(D == 1) {
    s.executor()
      .named("external_source")
      .forall(i, (m.template cells<ax::x, dm::quantities>())) {

      // Explicitly updating conserved variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i) += mass_density(i) * ga(i);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i) +=
        velocity(i).x() * mass_density(i) * ga(i).x();
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().named("external_source").forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      // Explicitly updating conserved variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j) += mass_density(i, j) * ga(i, j);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i, j) +=
        velocity(i, j).x() * ga(i, j).x() +
        velocity(i, j).y() * mass_density(i, j) * ga(i, j).y();
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().named("external_source").forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      // Explicitly updating conserved variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j, k) += mass_density(i, j, k) * ga(i, j, k);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i, j, k) +=
        velocity(i, j, k).x() * mass_density(i, j, k) * ga(i, j, k).x() +
        velocity(i, j, k).y() * mass_density(i, j, k) * ga(i, j, k).y() +
        velocity(i, j, k).z() * mass_density(i, j, k) * ga(i, j, k).z();
    }; // forall
  }
} // explicitSourceUpdate

} // namespace hard::tasks

#endif // HARD_MODULE_HYDRO_EXTERNAL_SOURCE_HH
