#ifndef HARD_MODULES_RAD_TASKS_INITIAL_KELVIN_HELM_RAD_HH
#define HARD_MODULES_RAD_TASKS_INITIAL_KELVIN_HELM_RAD_HH

#include "../../constants.hh"
#include "options.hh"
#include "spec/types.hh"
#include "types.hh"
#include <../modules/spec/utils.hh>

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
  spec::config_py config(opt::config.value());

  // Domain Parameters -
  // X and Z dimensions not used currently so not called here
  const double y_min = config["coords"][0][1].cast<double>();
  const double y_max = config["coords"][1][1].cast<double>();
  const double L_y = std::abs(y_max - y_min);

  // setting density and velocity
  // H is the heavier fluid at bottom
  const double rH =
    config["problem_parameters"]["fluid_mass_density_high"].cast<double>();
  const double pH =
    config["problem_parameters"]["fluid_pressure_high"].cast<double>();
  const double uH =
    config["problem_parameters"]["fluid_x_velocity_high"].cast<double>();
  const double vH =
    config["problem_parameters"]["fluid_y_velocity_high"].cast<double>();
  double w_h{0};

  // L is the lighter fluid at top
  const double rL =
    config["problem_parameters"]["fluid_mass_density_low"].cast<double>();
  const double pL =
    config["problem_parameters"]["fluid_pressure_low"].cast<double>();
  const double uL =
    config["problem_parameters"]["fluid_x_velocity_low"].cast<double>();
  const double vL =
    config["problem_parameters"]["fluid_y_velocity_low"].cast<double>();
  double w_l{0};

  if constexpr(D == 3) {
    // setting density and velocity
    // H is the heavier fluid at bottom
    w_h = config["problem_parameters"]["fluid_z_velocity_high"].cast<double>();
    // L is the lighter fluid at top
    w_l = config["problem_parameters"]["fluid_z_velocity_low"].cast<double>();
  }

  // setting fluid separation and velocity perturbation fractions and
  // location(s) (fraction so cast to make actual location agnostic of domain
  // dimensions)
  const double y_sep_f =
    config["problem_parameters"]["fluid_sep_y_frac"].cast<double>();
  const double y_vp_f =
    config["problem_parameters"]["vel_per_y_frac"].cast<double>();
  const double y_sep = y_sep_f * L_y;
  const double y_vp = y_vp_f * L_y;

  // setting perturbation parameters
  //  constants for perturbation
  const double N = config["problem_parameters"]["perturb_N"].cast<double>();

  const double wavenumber = 2.0 * M_PI * N;

  // setting temperature
  double rad_temp = config["problem_parameters"]["rad_temp"].cast<double>();

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented");
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
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
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);

          // initialize two different density and velocity fluids
          if(std::abs(y - (L_y / 2)) > y_sep) {
            mass_density(i, j, k) = rL;
            momentum_density(i, j, k).x() = rL * uL;
            momentum_density(i, j, k).y() = rL * vL;
            momentum_density(i, j, k).z() = rL * w_l;
            const double e = util::find_sie(eos, rL, pL);
            total_energy_density(i, j, k) = rL * e + 0.5 * rL * (vL * vL);
          }
          else {
            mass_density(i, j, k) = rH;
            momentum_density(i, j, k).x() = rH * uH;
            momentum_density(i, j, k).y() = rH * vH;
            momentum_density(i, j, k).z() = rH * w_h;
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

#endif // HARD_MODULES_RAD_TASKS_INITIAL_KELVIN_HELM_RAD_HH
