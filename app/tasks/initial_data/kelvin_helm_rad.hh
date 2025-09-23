#pragma once

#include "../../constants.hh"
#include "../../options.hh"
#include "../../types.hh"
#include <cmath>
#include <cstddef>
#include <spec/utils.hh>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

/*----------------------------------------------------------------------------*
  Kelvin-Helmholtz Instability with Radiation Enabled
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
kh_instability_rad(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  const eos::eos_wrapper & eos) noexcept {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // constants required from constants.hh
  const double a = hard::constants::cgs::radiation_constant;

  // Parse input parameters
  YAML::Node config = YAML::LoadFile(opt::config.value());

  // Domain Parameters -
  // X and Z dimensions not used currently so not called here
  const double y_min = config["coords"][0][1].as<double>();
  const double y_max = config["coords"][1][1].as<double>();
  const double L_y = std::abs(y_max - y_min);

  // setting density and velocity
  // H is the heavier fluid at bottom
  const double rH =
    config["problem_parameters"]["fluid_mass_density_high"].as<double>();
  const double pH =
    config["problem_parameters"]["fluid_pressure_high"].as<double>();
  const double uH =
    config["problem_parameters"]["fluid_x_velocity_high"].as<double>();
  const double vH =
    config["problem_parameters"]["fluid_y_velocity_high"].as<double>();
  double wH{0};

  // L is the lighter fluid at top
  const double rL =
    config["problem_parameters"]["fluid_mass_density_low"].as<double>();
  const double pL =
    config["problem_parameters"]["fluid_pressure_low"].as<double>();
  const double uL =
    config["problem_parameters"]["fluid_x_velocity_low"].as<double>();
  const double vL =
    config["problem_parameters"]["fluid_y_velocity_low"].as<double>();
  double wL{0};

  if constexpr(D == 3) {
    // setting density and velocity
    // H is the heavier fluid at bottom
    wH = config["problem_parameters"]["fluid_z_velocity_high"].as<double>();
    // L is the lighter fluid at top
    wL = config["problem_parameters"]["fluid_z_velocity_low"].as<double>();
  }

  // setting fluid separation and velocity perturbation fractions and
  // location(s) (fraction so as to make actual location agnostic of domain
  // dimensions)
  const double y_sep_f =
    config["problem_parameters"]["fluid_sep_y_frac"].as<double>();
  const double y_vp_f =
    config["problem_parameters"]["vel_per_y_frac"].as<double>();
  const double y_sep = y_sep_f * L_y;
  const double y_vp = y_vp_f * L_y;

  // setting perturbation parameters
  //  constants for perturbation
  const double N = config["problem_parameters"]["perturb_N"].as<double>();

  const double wavenumber = 2.0 * M_PI * N;

  // setting temperature
  double rad_temp = 0;
#ifdef ENABLE_RADIATION
  rad_temp = config["problem_parameters"]["rad_temp"].as<double>();
#endif

  if constexpr(D == 1) {
    assert(
      false &&
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented");
  }
  else if constexpr(D == 2) {
    s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        // initialize two different density and velocity fluids
        if(std::abs(y - (L_y / 2)) > y_sep) {
          mass_density(i, j) = rL;
          momentum_density(i, j).x() = rL * uL;
          momentum_density(i, j).y() = rL * vL;
          const double e = util::find_sie(eos, rL, pL);
          total_energy_density(i, j) = rL * e + 0.5 * rL * (vL * vL);
        }
        else {
          mass_density(i, j) = rH;
          momentum_density(i, j).x() = rH * uH;
          momentum_density(i, j).y() = rH * vH;
          const double e = util::find_sie(eos, rH, pH);
          total_energy_density(i, j) = rH * e + 0.5 * rH * (vH * vH);
        } // if

        radiation_energy_density(i, j) =
          a * spec::utils::sqr(spec::utils::sqr(rad_temp));

        // velocity perturbations in the Y-direction
        if(std::abs(y - (L_y / 4)) < y_vp) {
          momentum_density(i, j).y() = 0.05 * sin(wavenumber * x);
        }
        if(std::abs(y - (3 * L_y / 4)) < y_vp) {
          momentum_density(i, j).y() = 0.05 * sin(wavenumber * x);
        }
      } // for
    }; // forall
  }
  else /* D == 3 */ {
    s.executor().forall(k, (m.template cells<ax::z, dm::quantities>())) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);

          // initialize two different density and velocity fluids
          if(std::abs(y - (L_y / 2)) > y_sep) {
            mass_density(i, j, k) = rL;
            momentum_density(i, j, k).x() = rL * uL;
            momentum_density(i, j, k).y() = rL * vL;
            momentum_density(i, j, k).z() = rL * wL;
            const double e = util::find_sie(eos, rL, pL);
            total_energy_density(i, j, k) = rL * e + 0.5 * rL * (vL * vL);
          }
          else {
            mass_density(i, j, k) = rH;
            momentum_density(i, j, k).x() = rH * uH;
            momentum_density(i, j, k).y() = rH * vH;
            momentum_density(i, j, k).z() = rH * wH;
            const double e = util::find_sie(eos, rH, pH);
            total_energy_density(i, j, k) = rH * e + 0.5 * rH * (vH * vH);
          } // if

          radiation_energy_density(i, j, k) =
            a * spec::utils::sqr(spec::utils::sqr(rad_temp));

          // velocity perturbations in the Y-direction
          if(std::abs(y - (L_y / 4)) < y_vp) {
            momentum_density(i, j, k).y() = 0.05 * sin(wavenumber * x);
          }
          if(std::abs(y - (3 * L_y / 4)) < y_vp) {
            momentum_density(i, j, k).y() = 0.05 * sin(wavenumber * x);
          }
        } // for
      }
    }; // forall
  } // if
} //  kh_instability

} // namespace hard::tasks::initial_data
