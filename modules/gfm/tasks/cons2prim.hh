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
  typename field<vec<Dim>>::template accessor<ro, na> momentum_density_a,
  field<double>::accessor<ro, na> total_energy_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> velocity_a,
  field<double>::accessor<wo, na> pressure_a,
  field<double>::accessor<wo, na> specific_internal_energy_a,
  field<double>::accessor<wo, na> soundspeed_a,
  typename field<int>::accessor<ro, ro> material_a,
  std::vector<eos::eos_wrapper> const & eos_mm) noexcept {

  eos::eos_dispatch<flecsi::exec::accelerator> eos_disp(eos_mm);

  /**
   * Important: we are using a flat iterator over the entire domain
   * including ghost cells. At this stage ghost cells might be not
   * correctly populated, triggering wrong const auto & eos = eos_disp[mat]
   * accessing
   */
  s.executor().forall(
    i, flecsi::util::iota_view({}, pressure_a.span().size())) {
    const auto mat = static_cast<std::size_t>(material_a(i));
    if(mat < eos_disp.size_) {
      const auto & eos = eos_disp[mat];
      velocity_a(i) = momentum_density_a(i) / mass_density_a(i);
      specific_internal_energy_a(i) =
        (total_energy_density_a(i) -
          0.5 * mass_density_a(i) * (velocity_a(i).norm_squared())) /
        mass_density_a(i);
      pressure_a(i) =
        eos.pRhoSie(mass_density_a(i), specific_internal_energy_a(i));
      soundspeed_a(i) =
        eos.cRhoSie(mass_density_a(i), specific_internal_energy_a(i));
    }
  };
}

} // namespace hard::tasks::gfm

#endif // HARD_MODULE_HYDRO_CONS2PRIM_HH
