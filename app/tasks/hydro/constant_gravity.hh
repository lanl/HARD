
#pragma once

#include "../../types.hh"
#include <cstddef>

namespace flastro::tasks::hydro {

template<std::size_t D>
void
applyConstantGravity(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  typename field<vec<D>>::template accessor<ro, na> ru_a,
  typename field<vec<D>>::template accessor<rw, na> dt_momentum_density_a,
  field<double>::accessor<rw, na> dt_energy_density_a) {
  auto mass_density = m.template mdcolex<is::cells>(r_a);
  auto momentum_density = m.template mdcolex<is::cells>(ru_a);
  auto dt_mom = m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_e = m.template mdcolex<is::cells>(dt_energy_density_a);

  const double g = 1.0;

  if constexpr(D == 3) {

    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          dt_e(i, j, k) += -g * momentum_density(i, j, k).z;
          dt_mom(i, j, k).z += -g * mass_density(i, j, k);
        }
      }
    }
  }
  else {
    flog_fatal("")
  }
}

} // namespace flastro::tasks::hydro
