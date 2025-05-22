
#pragma once

#include "../../constants.hh"
#include "../../types.hh"
#include "spec/utils.hh"
#include <cmath>
#include <cstddef>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

//
// A sedov blast wave is being initialized.
// A hotspot of high temperature is added at the center of a
// domain with constant density, zero velocity, and low temperature.
//
template<std::size_t Dim>
auto
sedov_blast(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const double x0 =
    config["problem_parameters"]["hotspot_position"][0][0].as<double>();
  const double y0 =
    config["problem_parameters"]["hotspot_position"][0][1].as<double>();
  const double z0 =
    config["problem_parameters"]["hotspot_position"][0][2].as<double>();
  const double radius =
    config["problem_parameters"]["hotspot_radius"].as<double>();
  const double E_0 = config["problem_parameters"]["E_0"].as<double>();

  // some problem parameters: for now hardwired
  const double density = 1.0;

  const double volume =
    pow(radius, Dim) * (Dim > 1 ? M_PI : 1.) * (Dim > 2 ? 4. / 3. : 1.);
  const double rE_inside = E_0 / volume;
  const double rE_outside = 1.0e-5;

  using spec::utils::sqr;

  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "init_sedov_1d") {
      const auto x = m.template center<ax::x>(i);

      double distance = std::abs(x - x0);

      mass_density(i) = density;
      momentum_density(i).x = 0.0; // velocity is zero.
      radiation_energy_density(i) = 0.0;
      if(distance < radius) {
        total_energy_density(i) = rE_inside;
      }
      else {
        total_energy_density(i) = rE_outside;
      }
    }; // forall
  }
  else if constexpr(Dim == 2) {
    forall(j, (m.template cells<ax::y, dm::quantities>()), "init_sedov_2d") {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        double distance = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));

        mass_density(i, j) = density;
        momentum_density(i, j).x = 0.0;
        momentum_density(i, j).y = 0.0;
        radiation_energy_density(i, j) = 0.0;
        if(distance < radius) {
          total_energy_density(i, j) = rE_inside;
        }
        else {
          total_energy_density(i, j) = rE_outside;
        }
      } // for
    }; // forall
  }
  else if constexpr(Dim == 3) {
    forall(k, (m.template cells<ax::z, dm::quantities>()), "init_sedov_3d") {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);
          const auto z = m.template center<ax::z>(k);

          double distance = std::sqrt(
            (x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0));

          mass_density(i, j, k) = density;
          momentum_density(i, j, k).x = 0.0;
          momentum_density(i, j, k).y = 0.0;
          momentum_density(i, j, k).z = 0.0;
          radiation_energy_density(i, j, k) = 0.0;
          if(distance < radius) {
            total_energy_density(i, j, k) = rE_inside;
          }
          else {
            total_energy_density(i, j, k) = rE_outside;
          }
        } // for
      }
    }; // forall
  }
}

} // namespace hard::tasks::initial_data
