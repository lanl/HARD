#ifndef HARD_HYDRO_MODULES_SEDOV_HH
#define HARD_HYDRO_MODULES_SEDOV_HH

#include "../../constants.hh"
#include "spec/types.hh"
#include "spec/utils.hh"

namespace hard::tasks::initial_data {

//
// A sedov blast wave is being initialized.
// A hotspot of high temperature is added at the center of a
// domain with constant density, zero velocity, and low temperature.
//
template<std::size_t D>
auto
sedov_blast(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<D>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);

  // Parse input parameters
  spec::config_py config(opt::config.value());

  const double x0 =
    config["problem_parameters"]["hotspot_position"][0][0].cast<double>();
  const double y0 =
    config["problem_parameters"]["hotspot_position"][0][1].cast<double>();
  const double z0 =
    config["problem_parameters"]["hotspot_position"][0][2].cast<double>();
  const double radius =
    config["problem_parameters"]["hotspot_radius"].cast<double>();
  const double E_0 = config["problem_parameters"]["E_0"].cast<double>();

  // some problem parameters: for now hardwired
  const double density = 1.0;

  const double volume =
    pow(radius, D) * (D > 1 ? M_PI : 1.) * (D > 2 ? 4. / 3. : 1.);
  const double rE_inside = E_0 / volume;
  const double rE_outside = 1.0e-5;

  using spec::utils::sqr;

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const auto x = m.template center<ax::x>(i);
      double distance = std::abs(x - x0);
      mass_density(i) = density;
      momentum_density(i).x() = 0.0;
      total_energy_density(i) = distance < radius ? rE_inside : rE_outside;
    }; // forall
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);
        double distance = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
        mass_density(i, j) = density;
        momentum_density(i, j) = vec<D>(0.0);
        total_energy_density(i, j) = distance < radius ? rE_inside : rE_outside;
      } // for
    }; // forall
  }
  else if constexpr(D == 3) {
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);
          const auto z = m.template center<ax::z>(k);
          double distance = std::sqrt(
            (x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0));
          mass_density(i, j, k) = density;
          momentum_density(i, j, k) = vec<D>(0.0);
          total_energy_density(i, j, k) =
            distance < radius ? rE_inside : rE_outside;
        } // for
      }
    }; // forall
  }
}

} // namespace hard::tasks::initial_data

#endif // HARD_HYDRO_MODULES_SEDOV_HH
