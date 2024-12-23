
#pragma once

#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

//
// The Liska-Wendroff implosion test
// (http://www-troja.fjfi.cvut.cz/~liska/CompareEuler/compare8/)
//
// This test problem runs on a 2D square domain [0, 0.3]^2. In the region x + y
// < 0.15, density=0.125 and pressure=0.14. Otherwise denisty and pressure are
// both set to 1.0. Initial velocity is zero everywhere.
//
template<std::size_t Dim>
auto
lw_implosion(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a) {

  if constexpr(Dim == 2) {

    auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
    auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
    auto total_energy_density =
      m.template mdcolex<is::cells>(total_energy_density_a);
    auto radiation_energy_density =
      m.template mdcolex<is::cells>(radiation_energy_density_a);

    auto const gamma = *gamma_a;
    const double mult = 1.0 / (gamma - 1.0);

    forall(j, (m.template cells<ax::y, dm::quantities>()), "init_sedov_2d") {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        if(x + y > 0.15) {
          mass_density(i, j) = 1.0;
          total_energy_density(i, j) = mult * 1.0;
        }
        else {
          mass_density(i, j) = 0.125;
          total_energy_density(i, j) = mult * 0.14;
        }

        momentum_density(i, j).x = 0.0;
        momentum_density(i, j).y = 0.0;

        radiation_energy_density(i, j) = 0.0;

      } // for
    }; // forall
  }
  else {
    flog_fatal("Liska-Wendroff implosion is a 2D problem")
  }
}

} // namespace hard::tasks::initial_data
