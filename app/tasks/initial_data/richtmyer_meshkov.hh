#pragma once

#include "../../constants.hh"
#include "../../options.hh"
#include "../../types.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

/*----------------------------------------------------------------------------*
  Richtmyer-Meshkov Instability.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
richtmyer_meshkov(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  const eos::eos_wrapper & eos) {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  YAML::Node config = YAML::LoadFile(opt::config.value());

  const double rAs =
    config["problem_parameters"]["density_above_shock"].as<double>();
  const double pAs =
    config["problem_parameters"]["pressure_above_shock"].as<double>();
  const double vAs =
    config["problem_parameters"]["velocity_above_shock"].as<double>();
  const double rL = config["problem_parameters"]["density_low"].as<double>();
  const double rH = config["problem_parameters"]["density_high"].as<double>();
  const double p0 = config["problem_parameters"]["pressure"].as<double>();
  const double amp = config["problem_parameters"]["perturbation"].as<double>();
  const double interface =
    config["problem_parameters"]["interface"].as<double>();
  // for single mode perturbation, wavelength = domain width.
  const double wavelength = config["coords"][1][0].as<double>();
  const double wavenumber = 2.0 * M_PI;

  // for defining shock
  const double y_shock = interface + 0.01;

  // constant for the radiation energy
  const double a = hard::constants::cgs::radiation_constant;

  // temperature for radiation
  const double rad_temp =
    config["radiation_parameters"]["rad_temp"].as<double>();
  flog(info) << "radiation temperature is " << rad_temp << std::endl;

  if constexpr(D == 2) {
    s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const double x = m.template center<ax::x>(i);
        const double y = m.template center<ax::y>(j);

        // define the initial perturbation at the interface
        const double interface_y =
          interface + amp * (1 - std::cos(wavenumber * x / wavelength)) / 2;

        // place holder values for imposing initial conditions
        double rho, vy = 0, vx = 0.0, p;

        // lighter fluid above shock
        if(y >= y_shock) {
          rho = rAs;
          p = pAs;
          vy = vAs;
        }
        // lighter fluid below shock
        else if(y >= interface_y && y < y_shock) {
          rho = rL; // 0.125
          p = p0; // 0.1
          vy = 0.0;
        }
        // heavy fluid
        else {
          rho = rH;
          p = p0;
          vy = 0.0;
        }
        // Assign fields
        mass_density(i, j) = rho;
        momentum_density(i, j) = vec<D>(0.0);
        momentum_density(i, j).x() = rho * vx;
        momentum_density(i, j).y() = rho * vy;
        const double e = util::find_sie(eos, rho, p);
        const double kinetic = 0.5 * rho * vy * vy;
        total_energy_density(i, j) = rho * e + kinetic;
        radiation_energy_density(i, j) = 0;

#ifdef ENABLE_RADIATION
        radiation_energy_density(i, j) =
          a * spec::utils::sqr(spec::utils::sqr(rad_temp));
#endif
      }
    };
  }
  else {
    flog_fatal("Richtmyer-Meshkov problem is only implemented for D == 2");
  }
}

} // namespace hard::tasks::initial_data
