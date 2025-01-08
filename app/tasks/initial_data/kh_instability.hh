
#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>
#include <spec/utils.hh>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

//
// Initial data to demonstrate the Kelvin-Helmholtz instability.
//

template<std::size_t D>
auto
kh_instability(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a) {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto const gamma = *gamma_a;

  const double mult = 1.0 / (gamma - 1.0);

  using spec::utils::sqr;

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const double background_density =
    config["problem_parameters"]["background_density"].as<double>();
  const double background_velocity =
    config["problem_parameters"]["background_velocity"].as<double>();
  const double strip_density =
    config["problem_parameters"]["strip_density"].as<double>();
  const double strip_velocity =
    config["problem_parameters"]["strip_velocity"].as<double>();
  const double pressure = config["problem_parameters"]["pressure"].as<double>();

  const size_t perturbation_mode =
    config["problem_parameters"]["perturbation_mode"].as<size_t>();
  const double perturbation_amplitude =
    config["problem_parameters"]["perturbation_amplitude"].as<double>();
  const double perturbation_width =
    config["problem_parameters"]["perturbation_width"].as<double>();

  const double wavenumber = perturbation_mode * 2.0 * M_PI;
  const double gaussian_factor = 0.5 / sqr(perturbation_width);

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    forall(j, (m.template cells<ax::y, dm::quantities>()), "init_kh_2d") {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        const double vy = perturbation_amplitude * sin(wavenumber * x) *
                          (exp(-gaussian_factor * sqr(y - 0.25)) +
                            exp(-gaussian_factor * sqr(y - 0.75)));

        // initialize two different density and velocity fluids
        if(std::abs(y - 0.5) > 0.25) { // background
          mass_density(i, j) = background_density;
          momentum_density(i, j).x = background_density * background_velocity;
          momentum_density(i, j).y = background_density * vy;
          total_energy_density(i, j) =
            mult * pressure +
            0.5 * background_density * (sqr(background_velocity) + sqr(vy));
        }
        else { // strip
          mass_density(i, j) = strip_density;
          momentum_density(i, j).x = strip_density * strip_velocity;
          momentum_density(i, j).y = strip_density * vy;
          total_energy_density(i, j) =
            mult * pressure +
            0.5 * strip_density * (sqr(strip_velocity) + sqr(vy));
        } // if

        radiation_energy_density(i, j) = 0.0;

      } // for
    }; // forall
  }
  else /* D == 3 */ {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 3 is not implemented")
  } // if
} //  kh_instability

} // namespace hard::tasks::initial_data
