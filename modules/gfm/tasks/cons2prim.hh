#ifndef HARD_MODULE_GFM_CONS2PRIM_HH
#define HARD_MODULE_GFM_CONS2PRIM_HH

#include "../modules/common/tasks/utils.hh"
#include "../modules/spec/eos.hh"
#include <cstddef>

namespace hard::tasks::gfm {

struct conservative_to_primitive_mm_task : annotation_task_region {
  inline static const std::string name{"conservative_to_primitive_mm"};
};

template<std::size_t Dim>
void
conservative_to_primitive_mm(flecsi::exec::accelerator s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<ro, na> mass_density_a,
  typename vector_field_double<Dim>::template accessor<ro, na>
    momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  typename vector_field_double<Dim>::template accessor<wo, na> velocity_a,
  field<double>::accessor<wo, na> pressure_a,
  field<double>::accessor<wo, na> specific_internal_energy_a,
  field<double>::accessor<wo, na> soundspeed_a,
  typename field<int>::accessor<ro, ro> material_a,
  std::vector<eos::eos_wrapper> const & eos_mm) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  eos::eos_dispatch<flecsi::exec::accelerator> eos_disp(eos_mm);

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density =
    vector_field_double<Dim>::mdcolex(m, momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto velocity = vector_field_double<Dim>::mdcolex(m, velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
  auto specific_internal_energy =
    m.template mdcolex<is::cells>(specific_internal_energy_a);
  auto soundspeed = m.template mdcolex<is::cells>(soundspeed_a);
  auto material = m.template mdcolex<is::cells>(material_a);

  /**
   * Important: we iterate over all cells including ghost cells.
   * Ghost cells might not be correctly populated, so we check material
   * validity.
   */
  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::all>())) {
      const auto mat = static_cast<std::size_t>(material(i));
      if(mat < eos_disp.size_) {
        const auto & eos = eos_disp[mat];
        velocity(i) = momentum_density(i) / mass_density(i);
        specific_internal_energy(i) =
          (total_energy_density(i) -
            0.5 * mass_density(i) * (velocity(i).norm_squared())) /
          mass_density(i);
        pressure(i) = eos.pRhoSie(mass_density(i), specific_internal_energy(i));
        soundspeed(i) =
          eos.cRhoSie(mass_density(i), specific_internal_energy(i));
      }
    }; // forall
  }
  else if constexpr(Dim == 2) {
    auto mdpolicy_yx = get_mdiota_policy(velocity[0],
      m.template cells<ax::y, dm::all>(),
      m.template cells<ax::x, dm::all>());

    s.executor().forall(ji, mdpolicy_yx) {
      auto [j, i] = ji;
      const auto mat = static_cast<std::size_t>(material(i, j));
      if(mat < eos_disp.size_) {
        const auto & eos = eos_disp[mat];
        velocity(i, j) = momentum_density(i, j) / mass_density(i, j);
        specific_internal_energy(i, j) =
          total_energy_density(i, j) / mass_density(i, j) -
          0.5 * velocity(i, j).norm_squared();
        pressure(i, j) =
          eos.pRhoSie(mass_density(i, j), specific_internal_energy(i, j));
        soundspeed(i, j) =
          eos.cRhoSie(mass_density(i, j), specific_internal_energy(i, j));
      }
    }; // forall
  }
  else /* Dim == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(velocity[0],
      m.template cells<ax::z, dm::all>(),
      m.template cells<ax::y, dm::all>(),
      m.template cells<ax::x, dm::all>());

    s.executor().forall(kji, mdpolicy_zyx) {
      auto [k, j, i] = kji;
      const auto mat = static_cast<std::size_t>(material(i, j, k));
      if(mat < eos_disp.size_) {
        const auto & eos = eos_disp[mat];
        velocity(i, j, k) = momentum_density(i, j, k) / mass_density(i, j, k);
        specific_internal_energy(i, j, k) =
          total_energy_density(i, j, k) / mass_density(i, j, k) -
          0.5 * velocity(i, j, k).norm_squared();
        pressure(i, j, k) =
          eos.pRhoSie(mass_density(i, j, k), specific_internal_energy(i, j, k));
        soundspeed(i, j, k) =
          eos.cRhoSie(mass_density(i, j, k), specific_internal_energy(i, j, k));
      }
    }; // forall
  } // if constexpr
}

} // namespace hard::tasks::gfm

#endif // HARD_MODULE_HYDRO_CONS2PRIM_HH
