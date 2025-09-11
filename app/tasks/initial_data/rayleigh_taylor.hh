#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

/*
Rayleigh-Taylor Instability
Reference: https://www.astro.princeton.edu/~jstone/Athena/tests/rt/rt.html
Currently this file is a copy of input file for Kelvin-Helmholtz Instability*/

template<std::size_t D>
auto
rt_instability(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  typename field<vec<D>>::template accessor<wo, na> gravity_force_a,
  typename single<vec<D>>::template accessor<ro> gravity_acc_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  const eos::eos_wrapper & eos) noexcept {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto gravity_force = m.template mdcolex<is::cells>(gravity_force_a);
  auto g = gravity_acc_a.get();
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  // fixed values to be used

  const double p0 =
    2.5; // For hydrostatic equilibrium pressure, add g*rho*y to p0
  const double v0 = 0.05; // For perturbed velocity, multiply
                          // ((1 + cos (4*pi*x))/2)*((1 + cos (3*pi*y))/2)
  // density and x-velocity for two fluids
  // L is the bottom fluid
  const double rL = 1.0;
  const double uL = 0.0;

  // H is the upper fluid
  const double rH = 2.0;
  const double uH = 0.0;

  // y-velocity and pressure for two fluids
  // v and p will be implemented as spatially varying arrays
  // vL, vH, pL, pH will therefore not be implemented
  // If needed, kindly uncomment this section and ensure correct values
  /*
  const double vL = v0;
  const double pL = p0;
  const double vH = v0;
  const double pH = p0;
  */

  // initialize variables to be used at location considered
  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {

      const auto x = m.template center<ax::x>(i);
      double r, u, v, p;
      // set the density and x-velocity based on y location
      if(x >= 0.75) {
        // upper heavy fluid
        r = rH;
        u = 0.0;
        p = p0 + rL * g.x() * 0.75 +
            g.x() * r * (x - 0.75); // Hydrostatic Equilibrium pressure
      }
      else {
        // bottom lighter fluid
        r = rL;
        u = 0.0;
        p = p0 + rL * g.x() * x;
      } // if

      // Gaussian-localized x-velocity perturbation
      v = v0 * 1 / 2 * std::exp(-std::pow((x - 0.75) / 0.05, 2));

      // Initialize state variables at location
      mass_density(i) = r;
      gravity_force(i).x() = r * g.x();
      momentum_density(i).x() = r * u;
      const double e = util::find_sie(eos, r, p);
      total_energy_density(i) = r * e + 0.5 * r * (u * u + v * v);
      radiation_energy_density(i) = 0.0;
    }; // forall
  }
  else if constexpr(D == 2) {
    s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        double r, u, v, p;
        // set the density and x-velocity based on y location
        if(y > 0.75) {
          // upper heavy fluid
          r = rH;
          u = uH;
          p = p0 + rL * g.y() * 0.75 +
              (g.y() * r * (y - 0.75)); // Hydrostatic Equilibrium pressure
        }
        else {
          // bottom lighter fluid
          r = rL;
          u = uL;
          p = p0 + rL * g.y() * y;
        } // if

        // initially perturbed y-velocity
        // v = v0 * ((1 + cos(4 * M_PI * (x - 0.25))) / 2) *
        //     ((1 + cos(3 * M_PI * (y - 0.75))) / 2);

        // Gaussian-localized y-velocity perturbation
        v = v0 * (1 + cos(4 * M_PI * (x - 0.25))) / 2 *
            std::exp(-std::pow((y - 0.75) / 0.05, 2));

        // Initialize state variables at location
        mass_density(i, j) = r;
        gravity_force(i, j).x() = r * g.x();
        gravity_force(i, j).y() = r * g.y();
        momentum_density(i, j).x() = r * u;
        momentum_density(i, j).y() = r * v;
        const double e = util::find_sie(eos, r, p);
        total_energy_density(i, j) = r * e + 0.5 * r * (u * u + v * v);
        radiation_energy_density(i, j) = 0.0;

      } // for
    }; // forall
  }
  else /* D == 3 */ {
    s.executor().forall(k, (m.template cells<ax::z, dm::quantities>())) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          const auto x = m.template center<ax::x>(i);
          const auto y = m.template center<ax::y>(j);
          const auto z = m.template center<ax::z>(k);

          double r, u, v, w, p;
          // set the density and x-velocity based on y location
          if(y > 0.75) {
            // upper heavy fluid
            r = rH;
            u = uH;
            w = 0.0;
            p = p0 + rL * g.y() * 0.75 +
                (g.y() * r * (y - 0.75)); // Hydrostatic Equilibrium pressure
          }
          else {
            // bottom lighter fluid
            r = rL;
            u = uL;
            w = 0.0;
            p = p0 + rL * g.y() * y;
          } // if

          // initially perturbed y-velocity
          // v = v0 * ((1 + cos(4 * M_PI * (x - 0.25))) / 2) *
          //     ((1 + cos(3 * M_PI * (y - 0.75))) / 2);

          // Gaussian-localized y-velocity perturbation
          v = v0 *
              (1 + cos(4 * M_PI * (x - 0.75)) * (cos(4 * M_PI * (z - 0.75)))) /
              2 * std::exp(-std::pow((y - 0.75) / 0.05, 2));

          // Initialize state variables at location
          mass_density(i, j, k) = r;
          gravity_force(i, j, k).x() = r * g.x();
          gravity_force(i, j, k).y() = r * g.y();
          gravity_force(i, j, k).z() = r * g.z();
          momentum_density(i, j, k).x() = r * u;
          momentum_density(i, j, k).y() = r * v;
          momentum_density(i, j, k).z() = r * w;
          const double e = util::find_sie(eos, r, p);
          total_energy_density(i, j, k) =
            r * e + 0.5 * r * (u * u + v * v + w * w);
          radiation_energy_density(i, j, k) = 0.0;
        } // for
      } // for
    }; // forall
  } // if
} //  rt_instability

} // namespace hard::tasks::initial_data
