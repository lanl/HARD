
#pragma once

#include "../../constants.hh"
#include "../../numerical_algorithms/root_finder.hh"
#include "../../types.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>
#include <spec/utils.hh>
#include <yaml-cpp/yaml.h>

namespace flastro::tasks::initial_data {

template<std::size_t D>
auto
radiative_kh_instability(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a,
  single<double>::accessor<ro> particle_mass_a) {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  // auto const gamma = *gamma_a;
  // auto const particle_mass = *particle_mass_a;

  const double gamma = 5.0 / 3.0;
  const double particle_mass = 1.67262192369e-24;

  using flastro::tasks::util::get_mdiota_policy;
  using spec::utils::sqr;

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const double y0 = config["problem_parameters"]["y0"].as<double>();

  // density, velocity and pressure for two strips of fluids
  const double strip_thickness =
    config["problem_parameters"]["strip_thickness"].as<double>();
  const double strip_mass_density =
    config["problem_parameters"]["strip_mass_density"].as<double>();
  const double strip_temperature =
    config["problem_parameters"]["strip_temperature"].as<double>();
  const double background_mass_density =
    config["problem_parameters"]["background_mass_density"].as<double>();
  const double velocity_shear =
    config["problem_parameters"]["velocity_shear"].as<double>();

  // Assuming P = 1/3 E for radiation part
  const double strip_gas_pressure = strip_mass_density *
                                    constants::cgs::boltzmann_constant *
                                    strip_temperature / particle_mass;
  const double strip_internal_energy_density =
    strip_gas_pressure / (gamma - 1.0);
  const double strip_radiation_energy_density =
    constants::cgs::radiation_constant * sqr(sqr(strip_temperature));
  const double pressure =
    strip_gas_pressure + (strip_radiation_energy_density / 3.0);

  const double c4 = constants::cgs::radiation_constant / 3.0;
  const double c1 = background_mass_density *
                    constants::cgs::boltzmann_constant / particle_mass;
  const double c0 = -pressure;
  const double background_temperature =
    numerical_algorithms::root_finder::halleys_get_root(
      c4, c1, c0, 1000 * strip_temperature);

  if(background_temperature < 0.0 or std::isnan(background_temperature)) {
    flog_fatal("Invalid value of T_bg (radiative_kh_instability). Please check "
               "the input parameters")
  }

  const double background_gas_pressure = background_mass_density *
                                         constants::cgs::boltzmann_constant *
                                         background_temperature / particle_mass;
  const double background_internal_energy_density =
    background_gas_pressure / (gamma - 1.0);
  const double background_radiation_energy_density =
    constants::cgs::radiation_constant * sqr(sqr(background_temperature));

  flog(info) << " -- Radiative KH Instability Problem -- " << std::endl;
  flog(info) << "  * Strip      : T = " << strip_temperature
             << ", P_gas = " << strip_gas_pressure
             << ", P_rad = " << strip_radiation_energy_density / 3.0 << "\n"
             << "  * Background : T = " << background_temperature
             << ", P_gas = " << background_gas_pressure
             << ", P_rad = " << background_radiation_energy_density / 3.0
             << std::endl;
  flog(info) << "  dP = "
             << strip_gas_pressure + strip_radiation_energy_density / 3.0 -
                  background_gas_pressure -
                  background_radiation_energy_density / 3.0
             << std::endl;

  // Parameters related to the velocity perturbation
  const double velocity_perturbation_amplitude =
    config["problem_parameters"]["velocity_perturbation_amplitude"]
      .as<double>();
  const double wavelength =
    config["problem_parameters"]["wavelength"].as<double>();
  const double perturb_width =
    config["problem_parameters"]["perturb_width"].as<double>();

  const double N = 2;
  const double wavenumber = N * 2.0 * M_PI / wavelength;
  const double gaussian_factor = 0.5 / sqr(perturb_width);

  const auto vx = 0.5 * velocity_shear;

  if constexpr(D == 3) {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_zyx, "init_radiative_kh") {

      auto [k, j, i] = kji;
      const auto x = m.template center<ax::x>(i);
      const auto y = m.template center<ax::y>(j);

      const double vy =
        velocity_perturbation_amplitude * sin(wavenumber * x) *
        (exp(-gaussian_factor * sqr(y - (y0 + 0.5 * strip_thickness))) +
          exp(-gaussian_factor * sqr(y - (y0 - 0.5 * strip_thickness))));

      if(std::abs(y - y0) > 0.5 * strip_thickness) {
        // background
        mass_density(i, j, k) = background_mass_density;
        momentum_density(i, j, k).x =
          -background_mass_density * 0.5 * velocity_shear;
        momentum_density(i, j, k).y = background_mass_density * vy;
        total_energy_density(i, j, k) =
          background_internal_energy_density +
          0.5 * background_mass_density * (sqr(vx) + sqr(vy));
        radiation_energy_density(i, j, k) = background_radiation_energy_density;
      }
      else {
        // strip at the center
        mass_density(i, j, k) = strip_mass_density;
        momentum_density(i, j, k).x = strip_mass_density * 0.5 * velocity_shear;
        momentum_density(i, j, k).y = strip_mass_density * vy;
        total_energy_density(i, j, k) =
          strip_internal_energy_density +
          0.5 * strip_mass_density * (sqr(vx) + sqr(vy));
        radiation_energy_density(i, j, k) = strip_radiation_energy_density;
      }

      // Z direction is dummy axis
      momentum_density(i, j, k).z = 0.0;
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_yx, "init_radiative_kh_2d") {

      auto [j, i] = ji;
      const auto x = m.template center<ax::x>(i);
      const auto y = m.template center<ax::y>(j);

      const double vy =
        velocity_perturbation_amplitude * sin(wavenumber * x) *
        (exp(-gaussian_factor * sqr(y - (y0 + 0.5 * strip_thickness))) +
          exp(-gaussian_factor * sqr(y - (y0 - 0.5 * strip_thickness))));

      if(std::abs(y - y0) > 0.5 * strip_thickness) {
        // background
        mass_density(i, j) = background_mass_density;
        momentum_density(i, j).x =
          -background_mass_density * 0.5 * velocity_shear;
        momentum_density(i, j).y = background_mass_density * vy;
        total_energy_density(i, j) =
          background_internal_energy_density +
          0.5 * background_mass_density * (sqr(vx) + sqr(vy));
        radiation_energy_density(i, j) = background_radiation_energy_density;
      }
      else {
        // strip at the center
        mass_density(i, j) = strip_mass_density;
        momentum_density(i, j).x = strip_mass_density * 0.5 * velocity_shear;
        momentum_density(i, j).y = strip_mass_density * vy;
        total_energy_density(i, j) =
          strip_internal_energy_density +
          0.5 * strip_mass_density * (sqr(vx) + sqr(vy));
        radiation_energy_density(i, j) = strip_radiation_energy_density;
      }

      const double y1 = y0 - 0.5 * strip_thickness;
      const double y2 = y0 + 0.5 * strip_thickness;
      const double transition_width = 0.02 * strip_thickness;

      double f = 0;
      if(std::abs(y - y0) < 0.5 * strip_thickness) {
        f = 1.0;
      }
      if(std::abs(y - y1) < 0.5 * transition_width) {
        f = 0.5 * (1 + sin(M_PI * (y - y1) / transition_width));
      }
      if(std::abs(y - y2) < 0.5 * transition_width) {
        f = 0.5 * (1 - sin(M_PI * (y - y2) / transition_width));
      }

      mass_density(i, j) =
        f * strip_mass_density + (1.0 - f) * background_mass_density;

      const double c4 = constants::cgs::radiation_constant / 3.0;
      const double c1 =
        mass_density(i, j) * constants::cgs::boltzmann_constant / particle_mass;
      const double c0 = -pressure;
      const double temperature =
        numerical_algorithms::root_finder::halleys_get_root(
          c4, c1, c0, 100 * strip_temperature);

      const double gas_pressure = mass_density(i, j) *
                                  constants::cgs::boltzmann_constant *
                                  temperature / particle_mass;
      const double internal_energy_density = gas_pressure / (gamma - 1.0);
      const double radE =
        constants::cgs::radiation_constant * sqr(sqr(temperature));

      double vx;
      if(std::abs(y - y0) > 0.5 * strip_thickness) {
        vx = -0.5 * velocity_shear;
      }
      else {
        vx = 0.5 * velocity_shear;
      }

      momentum_density(i, j).x = mass_density(i, j) * vx;
      momentum_density(i, j).y = mass_density(i, j) * vy;
      total_energy_density(i, j) =
        internal_energy_density +
        0.5 * mass_density(i, j) * (sqr(vx) + sqr(vy));
      radiation_energy_density(i, j) = radE;
    };
  }
  else {
    flog_fatal("Radiative KH instability problem for D == 1 is not implemented")
  }
} //  kh_instability

} // namespace flastro::tasks::initial_data
