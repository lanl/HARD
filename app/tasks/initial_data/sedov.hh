
#pragma once

#include "../../constants.hh"
#include "../../types.hh"
#include <cmath>
#include <cstddef>
#include <yaml-cpp/yaml.h>

namespace flastro::tasks::initial_data {

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
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a,
  single<double>::accessor<ro> particle_mass_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // adiabatic index, assuming ideal gas EOS
  auto const gamma = *gamma_a;
  auto const particle_mass = *particle_mass_a;
  const double mult = 1.0 / (gamma - 1.0);

  // physical constants in cgs units
  const double kb = flastro::constants::cgs::boltzmann_constant;
  const double a = flastro::constants::cgs::radiation_constant;

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

  // some problem parameters: for now hardwired
  const double density = 1.0e-2;
  const double temperature_inside = 8.0e6;
  const double temperature_outside = 1.0e2;

  if constexpr(Dim == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "init_sedov_1d") {
      const auto x = m.template center<ax::x>(i);

      double distance = std::abs(x - x0);

      mass_density(i) = density;
      momentum_density(i).x = 0.0; // velocity is zero.
      if(distance < radius) {
        double pressure = kb * temperature_inside * density / (particle_mass);
        total_energy_density(i) = mult * pressure;
        radiation_energy_density(i) = a * temperature_inside *
                                      temperature_inside * temperature_inside *
                                      temperature_inside;
      }
      else {
        double pressure = kb * temperature_outside * density / (particle_mass);
        total_energy_density(i) = mult * pressure;
        radiation_energy_density(i) = a * temperature_outside *
                                      temperature_outside *
                                      temperature_outside * temperature_outside;
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
        if(distance < radius) {
          double pressure = kb * temperature_inside * density / (particle_mass);
          total_energy_density(i, j) = mult * pressure;
          radiation_energy_density(i, j) =
            a * temperature_inside * temperature_inside * temperature_inside *
            temperature_inside;
        }
        else {
          double pressure =
            kb * temperature_outside * density / (particle_mass);
          total_energy_density(i, j) = mult * pressure;
          radiation_energy_density(i, j) =
            a * temperature_outside * temperature_outside *
            temperature_outside * temperature_outside;
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
          if(distance < radius) {
            double pressure =
              kb * temperature_inside * density / (particle_mass);
            total_energy_density(i, j, k) = mult * pressure;
            radiation_energy_density(i, j, k) =
              a * temperature_inside * temperature_inside * temperature_inside *
              temperature_inside;
          }
          else {
            double pressure =
              kb * temperature_outside * density / (particle_mass);
            total_energy_density(i, j, k) = mult * pressure;
            radiation_energy_density(i, j, k) =
              a * temperature_outside * temperature_outside *
              temperature_outside * temperature_outside;
          }
        } // for
      }
    }; // forall
  }
}

} // namespace flastro::tasks::initial_data
