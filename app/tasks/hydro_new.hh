

#pragma once

#include "../state.hh"
#include "hydro/compute_interface_fluxes.hh"
#include "hydro/reconstruct.hh"
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>

namespace flastro::tasks::hydro {

template<std::size_t Dim, typename Limiter>
void
compute_fluxes_terms(typename mesh<Dim>::template accessor<ro> m,
  // conserved variables to be updated
  field<double>::accessor<rw, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<rw, ro> momentum_density_a,
  field<double>::accessor<rw, ro> total_energy_density_a,
  field<double>::accessor<rw, ro> radiation_energy_density_a,
  // primitive varibles
  typename field<vec<Dim>>::template accessor<ro, ro> velocity_a,
  field<double>::accessor<ro, ro> pressure_a,
  // reconstructed primitives on faces
  field<double>::accessor<wo, ro> rTail_a,
  field<double>::accessor<wo, ro> rHead_a,
  typename field<vec<Dim>>::template accessor<wo, ro> uTail_a,
  typename field<vec<Dim>>::template accessor<wo, ro> uHead_a,
  field<double>::accessor<wo, ro> pTail_a,
  field<double>::accessor<wo, ro> pHead_a,
  field<double>::accessor<wo, ro> EradTail_a,
  field<double>::accessor<wo, ro> EradHead_a,
  // conservatives on faces, computed from the reconstructed primitives
  typename field<vec<Dim>>::template accessor<wo, ro> ruTail_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ruHead_a,
  field<double>::accessor<wo, ro> rETail_a,
  field<double>::accessor<wo, ro> rEHead_a,
  // Riemann fluxes at cell interfaces
  field<double>::accessor<wo, ro> rF_a,
  typename field<vec<Dim>>::template accessor<wo, ro> ruF_a,
  field<double>::accessor<wo, ro> rEF_a,
  field<double>::accessor<wo, ro> EradF_a,
  // Time derivative
  typename field<std::tuple<double, vec<Dim>, double, double>>::
    template accessor<wo, na> dudt_fluxes_a,
  //
  single<double>::accessor<ro> gamma_a,
  double dt) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
  auto rTail = m.template mdcolex<is::cells>(rTail_a);
  auto rHead = m.template mdcolex<is::cells>(rHead_a);
  auto uTail = m.template mdcolex<is::cells>(uTail_a);
  auto uHead = m.template mdcolex<is::cells>(uHead_a);
  auto pTail = m.template mdcolex<is::cells>(pTail_a);
  auto pHead = m.template mdcolex<is::cells>(pHead_a);
  auto ruTail = m.template mdcolex<is::cells>(ruTail_a);
  auto ruHead = m.template mdcolex<is::cells>(ruHead_a);
  auto rETail = m.template mdcolex<is::cells>(rETail_a);
  auto rEHead = m.template mdcolex<is::cells>(rEHead_a);
  auto EradTail = m.template mdcolex<is::cells>(EradTail_a);
  auto EradHead = m.template mdcolex<is::cells>(EradHead_a);
  auto rF = m.template mdcolex<is::cells>(rF_a);
  auto ruF = m.template mdcolex<is::cells>(ruF_a);
  auto rEF = m.template mdcolex<is::cells>(rEF_a);
  auto EradF = m.template mdcolex<is::cells>(EradF_a);
  auto dudt_fluxes = m.template mdcolex<is::cells>(dudt_fluxes_a);
  auto const gamma = *gamma_a;

  auto const one_over_gamma_minus_1 = 1.0 / (gamma - 1.0);

  using spec::utils::sqr;

  // Compute (dt / dx^i)
  const auto dt_over_dx_i = [&m, &dt]() {
    if constexpr(Dim == 1) {
      return std::array<double, 1>{dt / m.template delta<ax::x>()};
    }
    else if constexpr(Dim == 2) {
      return std::array<double, 2>{
        dt / m.template delta<ax::x>(), dt / m.template delta<ax::y>()};
    }
    else {
      return std::array<double, 3>{dt / m.template delta<ax::x>(),
        dt / m.template delta<ax::y>(),
        dt / m.template delta<ax::z>()};
    }
  }();

  if constexpr(Dim == 1) {

    // clang-format off
    compute_interface_fluxes<Dim, ax::axis::x>(m,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a,
      rF_a, ruF_a, rEF_a, EradF_a,
      gamma_a);
    // clang-format on

    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      get<0>(dudt_fluxes(i)) = -dt_over_dx_i.at(0) * (rF(i + 1) - rF(i));
      get<1>(dudt_fluxes(i)) = -dt_over_dx_i.at(0) * (ruF(i + 1) - ruF(i));
      get<2>(dudt_fluxes(i)) = -dt_over_dx_i.at(0) * (rEF(i + 1) - rEF(i));
#ifndef FLASTRO_DISABLE_RADIATION
      get<3>(dudt_fluxes(i)) = -dt_over_dx_i.at(0) * (EradF(i + 1) - EradF(i));
#endif
    }
  }
  else if constexpr(Dim == 2) {

    // clang-format off
    reconstruct<Dim, ax::axis::x, Limiter>(m,
      mass_density_a, velocity_a, pressure_a, radiation_energy_density_a,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a, gamma_a);
    compute_interface_fluxes<Dim, ax::axis::x>(m,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a,
      rF_a, ruF_a, rEF_a, EradF_a,
      gamma_a);
    // clang-format on

    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        get<0>(dudt_fluxes(i, j)) =
          -dt_over_dx_i.at(0) * (rF(i + 1, j) - rF(i, j));
        get<1>(dudt_fluxes(i, j)) =
          -dt_over_dx_i.at(0) * (ruF(i + 1, j) - ruF(i, j));
        get<2>(dudt_fluxes(i, j)) =
          -dt_over_dx_i.at(0) * (rEF(i + 1, j) - rEF(i, j));
#ifndef FLASTRO_DISABLE_RADIATION
        get<3>(dudt_fluxes(i, j)) =
          -dt_over_dx_i.at(0) * (EradF(i + 1, j) - EradF(i, j));
#endif
      }
    }

    // clang-format off
    reconstruct<Dim, ax::axis::y, Limiter>(m,
      mass_density_a, velocity_a, pressure_a, radiation_energy_density_a,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a, gamma_a);
    compute_interface_fluxes<Dim, ax::axis::y>(m,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a,
      rF_a, ruF_a, rEF_a, EradF_a,
      gamma_a);
    // clang-format on

    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        get<0>(dudt_fluxes(i, j)) +=
          -dt_over_dx_i.at(1) * (rF(i, j + 1) - rF(i, j));
        get<1>(dudt_fluxes(i, j)) +=
          -dt_over_dx_i.at(1) * (ruF(i, j + 1) - ruF(i, j));
        get<2>(dudt_fluxes(i, j)) +=
          -dt_over_dx_i.at(1) * (rEF(i, j + 1) - rEF(i, j));
#ifndef FLASTRO_DISABLE_RADIATION
        get<3>(dudt_fluxes(i, j)) +=
          -dt_over_dx_i.at(1) * (EradF(i, j + 1) - EradF(i, j));
#endif
      }
    }
  }
  else { // Dim == 3

    // clang-format off
    reconstruct<Dim, ax::axis::x, Limiter>(m,
      mass_density_a, velocity_a, pressure_a, radiation_energy_density_a,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a, gamma_a);
    compute_interface_fluxes<Dim, ax::axis::x>(m,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a,
      rF_a, ruF_a, rEF_a, EradF_a,
      gamma_a);
    // clang-format on

    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          get<0>(dudt_fluxes(i, j, k)) =
            -dt_over_dx_i.at(0) * (rF(i + 1, j, k) - rF(i, j, k));
          get<1>(dudt_fluxes(i, j, k)) =
            -dt_over_dx_i.at(0) * (ruF(i + 1, j, k) - ruF(i, j, k));
          get<2>(dudt_fluxes(i, j, k)) =
            -dt_over_dx_i.at(0) * (rEF(i + 1, j, k) - rEF(i, j, k));
#ifndef FLASTRO_DISABLE_RADIATION
          get<3>(dudt_fluxes(i, j, k)) =
            -dt_over_dx_i.at(0) * (EradF(i + 1, j, k) - EradF(i, j, k));
#endif
        }
      }
    }

    // clang-format off
    reconstruct<Dim, ax::axis::y, Limiter>(m,
      mass_density_a, velocity_a, pressure_a, radiation_energy_density_a,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a, gamma_a);
    compute_interface_fluxes<Dim, ax::axis::y>(m,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a,
      rF_a, ruF_a, rEF_a, EradF_a,
      gamma_a);
    // clang-format on

    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          get<0>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(1) * (rF(i, j + 1, k) - rF(i, j, k));
          get<1>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(1) * (ruF(i, j + 1, k) - ruF(i, j, k));
          get<2>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(1) * (rEF(i, j + 1, k) - rEF(i, j, k));
#ifndef FLASTRO_DISABLE_RADIATION
          get<3>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(1) * (EradF(i, j + 1, k) - EradF(i, j, k));
#endif
        }
      }
    }

    // clang-format off
    reconstruct<Dim, ax::axis::z, Limiter>(m,
      mass_density_a, velocity_a, pressure_a, radiation_energy_density_a,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a, gamma_a);
    compute_interface_fluxes<Dim, ax::axis::z>(m,
      rTail_a, rHead_a, uTail_a, uHead_a, pTail_a, pHead_a, EradTail_a, EradHead_a,
      ruTail_a, ruHead_a, rETail_a, rEHead_a,
      rF_a, ruF_a, rEF_a, EradF_a,
      gamma_a);
    // clang-format on

    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          get<0>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(2) * (rF(i, j, k + 1) - rF(i, j, k));
          get<1>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(2) * (ruF(i, j, k + 1) - ruF(i, j, k));
          get<2>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(2) * (rEF(i, j, k + 1) - rEF(i, j, k));
#ifndef FLASTRO_DISABLE_RADIATION
          get<3>(dudt_fluxes(i, j, k)) +=
            -dt_over_dx_i.at(2) * (EradF(i, j, k + 1) - EradF(i, j, k));
#endif
        }
      }
    }
  }
}

} // namespace flastro::tasks::hydro
