#ifndef GFM_TASKS_COMPUTE_SOUND_SPEED_HH
#define GFM_TASKS_COMPUTE_SOUND_SPEED_HH

#include "../modules/spec/eos.hh"
#include "spec/mesh.hh"
#include "spec/types.hh"

namespace hard::tasks::gfm {

struct compute_sound_speed_from_eos_task : annotation_task_region {
  inline static const std::string name{"compute_sound_speed_from_eos"};
};

template<std::size_t Dim>
void
compute_sound_speed_from_eos(flecsi::exec::accelerator s,
  typename mesh<Dim>::template accessor<ro> m_a,
  field<double>::accessor<ro, na> mass_density_a,
  field<double>::accessor<ro, na> specific_internal_energy_a,
  typename field<int>::accessor<ro, ro> material_a,
  field<double>::accessor<wo, na> soundspeed_a,
  std::vector<eos::eos_wrapper> const & eos_mm) noexcept {

  eos::eos_dispatch<flecsi::exec::accelerator> eos_disp(eos_mm);

  /**
   * Important: we are using a flat iterator over the entire domain
   * including ghost cells. At this stage ghost cells might be not
   * correctly populated, triggering wrong const auto & eos = eos_disp[mat]
   * accessing
   */
  s.executor().forall(
    i, flecsi::util::iota_view({}, soundspeed_a.span().size())) {
    const auto mat = static_cast<std::size_t>(material_a(i));
    if(mat < eos_disp.size_) {
      const auto & eos = eos_disp[mat];
      soundspeed_a(i) =
        eos.cRhoSie(mass_density_a(i), specific_internal_energy_a(i));
    }
  };
}

} // namespace hard::tasks::gfm

#endif