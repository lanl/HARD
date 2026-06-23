#ifndef GFM_TASKS_COMPUTE_SOUND_SPEED_HH
#define GFM_TASKS_COMPUTE_SOUND_SPEED_HH

#include "../modules/spec/eos.hh"
#include "spec/mesh.hh"
#include "spec/types.hh"

namespace hard::tasks::gfm {

template<std::size_t Dim>
void
compute_sound_speed_from_eos(flecsi::exec::accelerator s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<ro, na> mass_density,
  // field<double>::accessor<ro, na> pressure_a,
  field<double>::accessor<ro, na> specific_internal_energy,
  typename field<int>::accessor<ro, ro> material,
  field<double>::accessor<wo, na> soundspeed,
  std::vector<eos::eos_wrapper> const & eos_mm) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  // std::cout << "---- speed_sound ----" << std::endl;

  if constexpr(Dim == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      // velocity(i) = momentum_density(i) / mass_density(i);
      // specific_internal_energy(i) =
      //   (total_energy_density(i) -
      //     0.5 * mass_density(i) * (velocity(i).norm_squared())) /
      //   mass_density(i);
      const auto mat = static_cast<std::size_t>(material(i));
      const auto & eos = eos_mm[mat];

      // pressure(i) = eos.pRhoSie(mass_density(i),
      // specific_internal_energy(i));
      soundspeed(i) = eos.cRhoSie(mass_density(i), specific_internal_energy(i));
    }; // for
  }
} // if
} // namespace hard::tasks::gfm
  // namespace hard::tasks::gfm

#endif