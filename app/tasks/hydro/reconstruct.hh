
#pragma once

#include "../../types.hh"
#include <cstddef>

namespace hard::tasks::hydro {

//
// Perform reconstruction of primitive variables on cell interfaces, calculate
// corresponding conservative variables on faces, and store them into `*Head`
// and `*Tail` variables.
//
//
//    |      U_i     |
//    |       *      |
//     ^HEAD(i)     ^Tail(i)
//    |              |
//   ^Tail(i-1)       ^Head(i+1)
//
//
template<std::size_t Dim, typename Limiter>
void
reconstruct(std::size_t reconstruction_axis,
  typename mesh<Dim>::template accessor<ro> m,
  // cell-centered primitive varibles
  field<double>::accessor<ro, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, ro> velocity_a,
  field<double>::accessor<ro, ro> pressure_a,
  field<double>::accessor<ro, ro>
#ifndef DISABLE_RADIATION
    radiation_energy_density_a
#endif
  ,
  // reconstructed primitives on faces
  field<double>::accessor<wo, na> rTail_a,
  field<double>::accessor<wo, na> rHead_a,
  typename field<vec<Dim>>::template accessor<wo, na> uTail_a,
  typename field<vec<Dim>>::template accessor<wo, na> uHead_a,
  field<double>::accessor<wo, na> pTail_a,
  field<double>::accessor<wo, na> pHead_a,
  field<double>::accessor<wo, na>
#ifndef DISABLE_RADIATION
    EradTail_a
#endif
  ,
  field<double>::accessor<wo, na>
#ifndef DISABLE_RADIATION
    EradHead_a
#endif
  ,
  // reconstructed conservatives on faces
  typename field<vec<Dim>>::template accessor<wo, na> ruTail_a,
  typename field<vec<Dim>>::template accessor<wo, na> ruHead_a,
  field<double>::accessor<wo, na> rETail_a,
  field<double>::accessor<wo, na> rEHead_a,
  //
  single<double>::accessor<ro> gamma_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
#ifndef DISABLE_RADIATION
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
#endif
  auto rTail = m.template mdcolex<is::cells>(rTail_a);
  auto rHead = m.template mdcolex<is::cells>(rHead_a);
  auto uTail = m.template mdcolex<is::cells>(uTail_a);
  auto uHead = m.template mdcolex<is::cells>(uHead_a);
  auto pTail = m.template mdcolex<is::cells>(pTail_a);
  auto pHead = m.template mdcolex<is::cells>(pHead_a);
#ifndef DISABLE_RADIATION
  auto EradTail = m.template mdcolex<is::cells>(EradTail_a);
  auto EradHead = m.template mdcolex<is::cells>(EradHead_a);
#endif
  auto ruTail = m.template mdcolex<is::cells>(ruTail_a);
  auto ruHead = m.template mdcolex<is::cells>(ruHead_a);
  auto rETail = m.template mdcolex<is::cells>(rETail_a);
  auto rEHead = m.template mdcolex<is::cells>(rEHead_a);

  using hard::tasks::util::get_mdiota_policy;
  using spec::utils::sqr;

  if constexpr(Dim == 1) {
    // if(reconstruction_axis != 0) {
    //   flog_fatal("1D executable can only reconstruct along X axis");
    // }

    forall(i, (m.template cells<ax::x, dm::predictor>()), "reconstruct_1d") {
      auto const gamma = *gamma_a;
      auto const one_over_gamma_minus_1 = 1.0 / (gamma - 1.0);

      const auto r_face = Limiter::reconstruct(mass_density(i - 2),
        mass_density(i - 1),
        mass_density(i),
        mass_density(i + 1),
        mass_density(i + 2));
      const auto v_face = Limiter::reconstruct(velocity(i - 2).x,
        velocity(i - 1).x,
        velocity(i).x,
        velocity(i + 1).x,
        velocity(i + 2).x);
      const auto p_face = Limiter::reconstruct(pressure(i - 2),
        pressure(i - 1),
        pressure(i),
        pressure(i + 1),
        pressure(i + 2));

      // NOTE: Avoid at(i) for GPU. at() does some bound check which doesn't
      // work on device. It should though, since at is constexpr.
      rHead(i) = r_face[0];
      rTail(i) = r_face[1];
      uHead(i).x = v_face[0];
      uTail(i).x = v_face[1];
      pHead(i) = p_face[0];
      pTail(i) = p_face[1];

      // Compute conservative variables
      ruHead(i) = rHead(i) * uHead(i);
      ruTail(i) = rTail(i) * uTail(i);
      rEHead(i) = one_over_gamma_minus_1 * pHead(i) +
                  0.5 * rHead(i) * uHead(i).norm_squared();
      rETail(i) = one_over_gamma_minus_1 * pTail(i) +
                  0.5 * rTail(i) * uTail(i).norm_squared();

#ifndef DISABLE_RADIATION
      const auto Erad_face =
        Limiter::reconstruct(radiation_energy_density(i - 2),
          radiation_energy_density(i - 1),
          radiation_energy_density(i),
          radiation_energy_density(i + 1),
          radiation_energy_density(i + 2));
      EradHead(i) = Erad_face[0];
      EradTail(i) = Erad_face[1];
#endif
    }; // forall
  }
  else if constexpr(Dim == 2) {

    auto mdpolicy_pp = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    forall(ji, mdpolicy_pp, "reconstruct_2d") {
      auto const gamma = *gamma_a;
      auto const one_over_gamma_minus_1 = 1.0 / (gamma - 1.0);

      auto [j, i] = ji;
      std::array<double, 2> r_face;
      std::array<double, 2> vx_face;
      std::array<double, 2> vy_face;
      std::array<double, 2> p_face;

#ifndef DISABLE_RADIATION
      std::array<double, 2> Erad_face;
#endif

      if(reconstruction_axis == 0) {
        r_face = Limiter::reconstruct(mass_density(i - 2, j),
          mass_density(i - 1, j),
          mass_density(i, j),
          mass_density(i + 1, j),
          mass_density(i + 2, j));
        vx_face = Limiter::reconstruct(velocity(i - 2, j).x,
          velocity(i - 1, j).x,
          velocity(i, j).x,
          velocity(i + 1, j).x,
          velocity(i + 2, j).x);
        vy_face = Limiter::reconstruct(velocity(i - 2, j).y,
          velocity(i - 1, j).y,
          velocity(i, j).y,
          velocity(i + 1, j).y,
          velocity(i + 2, j).y);
        p_face = Limiter::reconstruct(pressure(i - 2, j),
          pressure(i - 1, j),
          pressure(i, j),
          pressure(i + 1, j),
          pressure(i + 2, j));

#ifndef DISABLE_RADIATION
        Erad_face = Limiter::reconstruct(radiation_energy_density(i - 2, j),
          radiation_energy_density(i - 1, j),
          radiation_energy_density(i, j),
          radiation_energy_density(i + 1, j),
          radiation_energy_density(i + 2, j));
#endif
      }
      else if(reconstruction_axis == 1) {
        r_face = Limiter::reconstruct(mass_density(i, j - 2),
          mass_density(i, j - 1),
          mass_density(i, j),
          mass_density(i, j + 1),
          mass_density(i, j + 2));
        vx_face = Limiter::reconstruct(velocity(i, j - 2).x,
          velocity(i, j - 1).x,
          velocity(i, j).x,
          velocity(i, j + 1).x,
          velocity(i, j + 2).x);
        vy_face = Limiter::reconstruct(velocity(i, j - 2).y,
          velocity(i, j - 1).y,
          velocity(i, j).y,
          velocity(i, j + 1).y,
          velocity(i, j + 2).y);
        p_face = Limiter::reconstruct(pressure(i, j - 2),
          pressure(i, j - 1),
          pressure(i, j),
          pressure(i, j + 1),
          pressure(i, j + 2));

#ifndef DISABLE_RADIATION
        Erad_face = Limiter::reconstruct(radiation_energy_density(i, j - 2),
          radiation_energy_density(i, j - 1),
          radiation_energy_density(i, j),
          radiation_energy_density(i, j + 1),
          radiation_energy_density(i, j + 2));
#endif
      }
      // else {
      //   flog_fatal("2D executable cannot reconstruct along Z axis");
      // }

      rHead(i, j) = r_face[0];
      rTail(i, j) = r_face[1];
      uHead(i, j).x = vx_face[0];
      uTail(i, j).x = vx_face[1];
      uHead(i, j).y = vy_face[0];
      uTail(i, j).y = vy_face[1];
      pHead(i, j) = p_face[0];
      pTail(i, j) = p_face[1];

      // Compute conservative variables
      ruHead(i, j) = rHead(i, j) * uHead(i, j);
      ruTail(i, j) = rTail(i, j) * uTail(i, j);
      rEHead(i, j) = one_over_gamma_minus_1 * pHead(i, j) +
                     0.5 * rHead(i, j) * uHead(i, j).norm_squared();
      rETail(i, j) = one_over_gamma_minus_1 * pTail(i, j) +
                     0.5 * rTail(i, j) * uTail(i, j).norm_squared();

#ifndef DISABLE_RADIATION
      EradHead(i, j) = Erad_face[0];
      EradTail(i, j) = Erad_face[1];
#endif
    }; // forall
  }
  else { // Dim == 3
    auto mdpolicy_ppp = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::predictor>(),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    forall(kji, mdpolicy_ppp, "reconstruct_3d") {
      auto const gamma = *gamma_a;
      auto const one_over_gamma_minus_1 = 1.0 / (gamma - 1.0);
      auto [k, j, i] = kji;
      std::array<double, 2> r_face;
      std::array<double, 2> vx_face;
      std::array<double, 2> vy_face;
      std::array<double, 2> vz_face;
      std::array<double, 2> p_face;

#ifndef DISABLE_RADIATION
      std::array<double, 2> Erad_face;
#endif

      if(reconstruction_axis == 0) {
        r_face = Limiter::reconstruct(mass_density(i - 2, j, k),
          mass_density(i - 1, j, k),
          mass_density(i, j, k),
          mass_density(i + 1, j, k),
          mass_density(i + 2, j, k));
        vx_face = Limiter::reconstruct(velocity(i - 2, j, k).x,
          velocity(i - 1, j, k).x,
          velocity(i, j, k).x,
          velocity(i + 1, j, k).x,
          velocity(i + 2, j, k).x);
        vy_face = Limiter::reconstruct(velocity(i - 2, j, k).y,
          velocity(i - 1, j, k).y,
          velocity(i, j, k).y,
          velocity(i + 1, j, k).y,
          velocity(i + 2, j, k).y);
        vz_face = Limiter::reconstruct(velocity(i - 2, j, k).z,
          velocity(i - 1, j, k).z,
          velocity(i, j, k).z,
          velocity(i + 1, j, k).z,
          velocity(i + 2, j, k).z);
        p_face = Limiter::reconstruct(pressure(i - 2, j, k),
          pressure(i - 1, j, k),
          pressure(i, j, k),
          pressure(i + 1, j, k),
          pressure(i + 2, j, k));

#ifndef DISABLE_RADIATION
        Erad_face = Limiter::reconstruct(radiation_energy_density(i - 2, j, k),
          radiation_energy_density(i - 1, j, k),
          radiation_energy_density(i, j, k),
          radiation_energy_density(i + 1, j, k),
          radiation_energy_density(i + 2, j, k));
#endif
      }
      else if(reconstruction_axis == 1) {
        r_face = Limiter::reconstruct(mass_density(i, j - 2, k),
          mass_density(i, j - 1, k),
          mass_density(i, j, k),
          mass_density(i, j + 1, k),
          mass_density(i, j + 2, k));
        vx_face = Limiter::reconstruct(velocity(i, j - 2, k).x,
          velocity(i, j - 1, k).x,
          velocity(i, j, k).x,
          velocity(i, j + 1, k).x,
          velocity(i, j + 2, k).x);
        vy_face = Limiter::reconstruct(velocity(i, j - 2, k).y,
          velocity(i, j - 1, k).y,
          velocity(i, j, k).y,
          velocity(i, j + 1, k).y,
          velocity(i, j + 2, k).y);
        vz_face = Limiter::reconstruct(velocity(i, j - 2, k).z,
          velocity(i, j - 1, k).z,
          velocity(i, j, k).z,
          velocity(i, j + 1, k).z,
          velocity(i, j + 2, k).z);
        p_face = Limiter::reconstruct(pressure(i, j - 2, k),
          pressure(i, j - 1, k),
          pressure(i, j, k),
          pressure(i, j + 1, k),
          pressure(i, j + 2, k));

#ifndef DISABLE_RADIATION
        Erad_face = Limiter::reconstruct(radiation_energy_density(i, j - 2, k),
          radiation_energy_density(i, j - 1, k),
          radiation_energy_density(i, j, k),
          radiation_energy_density(i, j + 1, k),
          radiation_energy_density(i, j + 2, k));
#endif
      }
      else { // Z-axis
        r_face = Limiter::reconstruct(mass_density(i, j, k - 2),
          mass_density(i, j, k - 1),
          mass_density(i, j, k),
          mass_density(i, j, k + 1),
          mass_density(i, j, k + 2));
        vx_face = Limiter::reconstruct(velocity(i, j, k - 2).x,
          velocity(i, j, k - 1).x,
          velocity(i, j, k).x,
          velocity(i, j, k + 1).x,
          velocity(i, j, k + 2).x);
        vy_face = Limiter::reconstruct(velocity(i, j, k - 2).y,
          velocity(i, j, k - 1).y,
          velocity(i, j, k).y,
          velocity(i, j, k + 1).y,
          velocity(i, j, k + 2).y);
        vz_face = Limiter::reconstruct(velocity(i, j, k - 2).z,
          velocity(i, j, k - 1).z,
          velocity(i, j, k).z,
          velocity(i, j, k + 1).z,
          velocity(i, j, k + 2).z);
        p_face = Limiter::reconstruct(pressure(i, j, k - 2),
          pressure(i, j, k - 1),
          pressure(i, j, k),
          pressure(i, j, k + 1),
          pressure(i, j, k + 2));

#ifndef DISABLE_RADIATION
        Erad_face = Limiter::reconstruct(radiation_energy_density(i, j, k - 2),
          radiation_energy_density(i, j, k - 1),
          radiation_energy_density(i, j, k),
          radiation_energy_density(i, j, k + 1),
          radiation_energy_density(i, j, k + 2));
#endif
      }

      rHead(i, j, k) = r_face[0];
      rTail(i, j, k) = r_face[1];
      uHead(i, j, k).x = vx_face[0];
      uTail(i, j, k).x = vx_face[1];
      uHead(i, j, k).y = vy_face[0];
      uTail(i, j, k).y = vy_face[1];
      uHead(i, j, k).z = vz_face[0];
      uTail(i, j, k).z = vz_face[1];
      pHead(i, j, k) = p_face[0];
      pTail(i, j, k) = p_face[1];

#ifndef DISABLE_RADIATION
      EradHead(i, j, k) = Erad_face[0];
      EradTail(i, j, k) = Erad_face[1];
#endif

      // Compute conservative variables on faces
      ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);
      ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
      rEHead(i, j, k) = one_over_gamma_minus_1 * pHead(i, j, k) +
                        0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();
      rETail(i, j, k) = one_over_gamma_minus_1 * pTail(i, j, k) +
                        0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();

    }; // forall
  }
}

} // namespace hard::tasks::hydro
