
#pragma once

#include "../../options.hh"
#include "../../types.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>
#include <cstring>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

//
// An acoustic wave set-up
//
template<std::size_t Dim>
auto
acoustic_wave(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  const eos::eos_wrapper & eos) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  YAML::Node config = YAML::LoadFile(opt::config.value());

  // Problem parameters
  // Equilibrium values
  const double r0{
    config["problem_parameters"]["r0"].as<double>()}; // Equilibrium density
  const double p0{
    config["problem_parameters"]["p0"].as<double>()}; // Equilibrium pressure

  // Perturbation amplitudes
  const double rA{
    config["problem_parameters"]["amplitude"].as<double>()}; // Density
  const double uA{
    config["problem_parameters"]["amplitude"].as<double>()}; // Velocity

  // Sound speed
  const double cs{sqrt(config["gamma"].as<double>() * p0 / r0)};

  //
  // Only 1D and 2D versions have been implemented.
  //
  if constexpr(Dim == 1) {

    // Define the wave number
    const double k{
      2 * M_PI * config["problem_parameters"]["scale"][0].as<double>()};

    // sine_quad is the volume average of the sine per cell, for sine wave
    // fvm initialization
    auto sine_quad = [k](double x0, double x1) {
      return (cos(k * x0) - cos(k * x1)) / (k * (x1 - x0));
    };

    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const auto x0{m.template head<ax::x>(i)};
      const auto x1{m.template tail<ax::x>(i)};
      const double ux{cs * uA * sine_quad(x0, x1)};
      mass_density(i) = r0 * (1 + rA * sine_quad(x0, x1));

      momentum_density(i).x() = mass_density(i) * ux;
      const double e = util::find_sie(eos, mass_density(i), p0);
      total_energy_density(i) =
        mass_density(i) * e + 0.5 * mass_density(i) * utils::sqr(ux);

      radiation_energy_density(i) = 0.0;

    }; // forall
  }
  else if constexpr(Dim == 2) {
    // Define the wave number
    const double kx{
      2 * M_PI * config["problem_parameters"]["scale"][0].as<double>()};
    const double ky{
      2 * M_PI * config["problem_parameters"]["scale"][1].as<double>()};
    const double k{std::sqrt(utils::sqr(kx) + utils::sqr(ky))};

    // sine_quad is the volume average of the sine per cell, for sine wave
    // fvm initialization
    auto sine_quad = [kx, ky](double x0, double x1, double y0, double y1) {
      return -((cos(kx * x1) - cos(kx * x0)) * (sin(ky * y1) - sin(ky * y0)) +
               (sin(kx * x1) - sin(kx * x0)) * (cos(ky * y1) - cos(ky * y0))) /
             (kx * ky * (x1 - x0) * (y1 - y0));
    };

    s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x0{m.template head<ax::x>(i)};
        const auto x1{m.template tail<ax::x>(i)};
        const auto y0{m.template head<ax::y>(j)};
        const auto y1{m.template tail<ax::y>(j)};

        const double ux{cs * uA * kx * sine_quad(x0, x1, y0, y1) / k};
        const double uy{cs * uA * ky * sine_quad(x0, x1, y0, y1) / k};
        mass_density(i, j) = r0 * (1 + rA * sine_quad(x0, x1, y0, y1));

        momentum_density(i, j).x() = mass_density(i, j) * ux;
        momentum_density(i, j).y() = mass_density(i, j) * uy;
        const double e = util::find_sie(eos, mass_density(i, j), p0);
        total_energy_density(i, j) =
          mass_density(i, j) * e +
          0.5 * mass_density(i, j) * (utils::sqr(ux) + utils::sqr(uy));

        radiation_energy_density(i, j) = 0.0;

      } // for
    }; // forall
  }
  else /* Dim == 3 */ {

    // Define the wave number
    const double kx{
      2 * M_PI * config["problem_parameters"]["scale"][0].as<double>()};
    const double ky{
      2 * M_PI * config["problem_parameters"]["scale"][1].as<double>()};
    const double kz{
      2 * M_PI * config["problem_parameters"]["scale"][2].as<double>()};
    const double k{std::sqrt(utils::sqr(kx) + utils::sqr(ky) + utils::sqr(kz))};

    // sine_quad is the volume average of the sine per cell, for sine wave
    // fvm initialization
    auto sine_quad =
      [kx, ky, kz](
        double x0, double x1, double y0, double y1, double z0, double z1) {
        return ((cos(kx * x1) - cos(kx * x0)) * (cos(ky * y1) - cos(ky * y0)) *
                   (cos(kz * z1) - cos(kz * z0)) -
                 (cos(kx * x1) - cos(kx * x0)) * (sin(ky * y1) - sin(ky * y0)) *
                   (sin(kz * z1) - sin(kz * z0)) -
                 (sin(kx * x1) - sin(kx * x0)) * (cos(ky * y1) - cos(ky * y0)) *
                   (sin(kz * z1) - sin(kz * z0)) -
                 (sin(kx * x1) - sin(kx * x0)) * (sin(ky * y1) - sin(ky * y0)) *
                   (cos(kz * z1) - cos(kz * z0))) /
               (kx * ky * kz * (x1 - x0) * (y1 - y0) * (z1 - z0));
      };

    s.executor().forall(l, (m.template cells<ax::z, dm::quantities>())) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          const auto x0{m.template head<ax::x>(i)};
          const auto x1{m.template tail<ax::x>(i)};
          const auto y0{m.template head<ax::y>(j)};
          const auto y1{m.template tail<ax::y>(j)};
          const auto z0{m.template head<ax::z>(l)};
          const auto z1{m.template tail<ax::z>(l)};

          const double ux{cs * uA * kx * sine_quad(x0, x1, y0, y1, z0, z1) / k};
          const double uy{cs * uA * ky * sine_quad(x0, x1, y0, y1, z0, z1) / k};
          const double uz{cs * uA * kz * sine_quad(x0, x1, y0, y1, z0, z1) / k};
          mass_density(i, j, l) =
            r0 * (1 + rA * sine_quad(x0, x1, y0, y1, z0, z1));

          momentum_density(i, j, l).x() = mass_density(i, j, l) * ux;
          momentum_density(i, j, l).y() = mass_density(i, j, l) * uy;
          momentum_density(i, j, l).z() = mass_density(i, j, l) * uz;
          const double e = util::find_sie(eos, mass_density(i, j, l), p0);
          total_energy_density(i, j, l) =
            mass_density(i, j, l) * e +
            0.5 * mass_density(i, j, l) *
              (utils::sqr(ux) + utils::sqr(uy) + utils::sqr(uz));

          radiation_energy_density(i, j, l) = 0.0;

        } // for
      } // for
    }; // forall
  }
} // acoustic_wave

} // namespace hard::tasks::initial_data
