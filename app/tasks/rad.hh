#ifndef HARD_TASKS_RAD_HH
#define HARD_TASKS_RAD_HH

#include "../constants.hh"
#include "../state.hh"
#include "utils.hh"
#include <cmath>
#include <cstddef>
#include <spec/utils.hh>

namespace hard::task::rad {

using hard::tasks::util::get_mdiota_policy;

// Get the gradient of the velocity using a 5-point stencil.
template<std::size_t D>
void
getGradV(typename mesh<D>::template accessor<ro> m,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    na> gradV_a,
  typename field<vec<D>>::template accessor<ro, ro> u_a) {

  auto u = m.template mdcolex<is::cells>(u_a);
  auto gradV = m.template mdcolex<is::cells>(gradV_a);

  if constexpr(D == 1) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());

    forall(i, (m.template cells<ax::x, dm::quantities>()), "getGradV") {
      gradV(i).xx =
        (u(i - 2).x - 8.0 * u(i - 1).x + 8.0 * u(i + 1).x - u(i + 2).x) *
        one_over_12dx;
    };
  }
  else if constexpr(D == 2) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());
    const double one_over_12dy = 1.0 / (12.0 * m.template delta<ax::y>());

    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "getGradV") {
      auto [j, i] = ji;

      gradV(i, j).xx = (u(i - 2, j).x - 8.0 * u(i - 1, j).x +
                         8.0 * u(i + 1, j).x - u(i + 2, j).x) *
                       one_over_12dx;
      gradV(i, j).xy = (u(i, j - 2).x - 8.0 * u(i, j - 1).x +
                         8.0 * u(i, j + 1).x - u(i, j + 2).x) *
                       one_over_12dy;

      gradV(i, j).yx = (u(i - 2, j).y - 8.0 * u(i - 1, j).y +
                         8.0 * u(i + 1, j).y - u(i + 2, j).y) *
                       one_over_12dx;
      gradV(i, j).yy = (u(i, j - 2).y - 8.0 * u(i, j - 1).y +
                         8.0 * u(i, j + 1).y - u(i, j + 2).y) *
                       one_over_12dy;
    };
  }
  else {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());
    const double one_over_12dy = 1.0 / (12.0 * m.template delta<ax::y>());
    const double one_over_12dz = 1.0 / (12.0 * m.template delta<ax::z>());

    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "getGradV") {
      auto [k, j, i] = kji;
      gradV(i, j, k).xx = (u(i - 2, j, k).x - 8.0 * u(i - 1, j, k).x +
                            8.0 * u(i + 1, j, k).x - u(i + 2, j, k).x) *
                          one_over_12dx;
      gradV(i, j, k).xy = (u(i, j - 2, k).x - 8.0 * u(i, j - 1, k).x +
                            8.0 * u(i, j + 1, k).x - u(i, j + 2, k).x) *
                          one_over_12dy;
      gradV(i, j, k).xz = (u(i, j, k - 2).x - 8.0 * u(i, j, k - 1).x +
                            8.0 * u(i, j, k + 1).x - u(i, j, k + 2).x) *
                          one_over_12dz;

      gradV(i, j, k).yx = (u(i - 2, j, k).y - 8.0 * u(i - 1, j, k).y +
                            8.0 * u(i + 1, j, k).y - u(i + 2, j, k).y) *
                          one_over_12dx;
      gradV(i, j, k).yy = (u(i, j - 2, k).y - 8.0 * u(i, j - 1, k).y +
                            8.0 * u(i, j + 1, k).y - u(i, j + 2, k).y) *
                          one_over_12dy;
      gradV(i, j, k).yz = (u(i, j, k - 2).y - 8.0 * u(i, j, k - 1).y +
                            8.0 * u(i, j, k + 1).y - u(i, j, k + 2).y) *
                          one_over_12dz;

      gradV(i, j, k).zx = (u(i - 2, j, k).z - 8.0 * u(i - 1, j, k).z +
                            8.0 * u(i + 1, j, k).z - u(i + 2, j, k).z) *
                          one_over_12dx;
      gradV(i, j, k).zy = (u(i, j - 2, k).z - 8.0 * u(i, j - 1, k).z +
                            8.0 * u(i, j + 1, k).z - u(i, j + 2, k).z) *
                          one_over_12dy;
      gradV(i, j, k).zz = (u(i, j, k - 2).z - 8.0 * u(i, j, k - 1).z +
                            8.0 * u(i, j, k + 1).z - u(i, j, k + 2).z) *
                          one_over_12dz;
    };
  }
} // getGradV

// Get the radiation pressure tensor P
template<std::size_t D>
void
getTensorP(typename mesh<D>::template accessor<ro> m,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    na> P_tensor_a,
  field<double>::accessor<ro, na> Esf_a,
  typename field<vec<D>>::template accessor<ro, na> gradEsf_a,
  field<double>::accessor<ro, na> gradE_mag_a,
  field<double>::accessor<ro, na> lambda_a,
  field<double>::accessor<ro, na> R_a) {

  auto P_tensor = m.template mdcolex<is::cells>(P_tensor_a);
  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto gradEsf = m.template mdcolex<is::cells>(gradEsf_a);
  auto gradE_mag = m.template mdcolex<is::cells>(gradE_mag_a);
  auto lambda = m.template mdcolex<is::cells>(lambda_a);
  auto R = m.template mdcolex<is::cells>(R_a);

  const double eps = 1.0e-15;

  const auto compute_eddington_factor = [](
                                          double lambda_value, double R_value) {
    return lambda_value + spec::utils::sqr(lambda_value * R_value);
  };

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "getTensorP") {
      const double f = compute_eddington_factor(lambda(i), R(i));
      P_tensor(i).xx = (0.5 * (1 - f) + 0.5 * (3 * f - 1)) * Esf(i);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(Esf,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "getTensorP") {
      auto [j, i] = ji;

      const double nx = gradEsf(i, j).x / (gradE_mag(i, j) + eps);
      const double ny = gradEsf(i, j).y / (gradE_mag(i, j) + eps);

      const double f = compute_eddington_factor(lambda(i, j), R(i, j));

      P_tensor(i, j).xx =
        (0.5 * (1 - f) + 0.5 * (3 * f - 1) * nx * nx) * Esf(i, j);
      P_tensor(i, j).xy = (0.5 * (3 * f - 1) * nx * ny) * Esf(i, j);
      P_tensor(i, j).yx = P_tensor(i, j).xy;
      P_tensor(i, j).yy =
        (0.5 * (1 - f) + 0.5 * (3 * f - 1) * ny * ny) * Esf(i, j);
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(Esf,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "getTensorP") {
      auto [k, j, i] = kji;

      const double nx = gradEsf(i, j, k).x / (gradE_mag(i, j, k) + eps);
      const double ny = gradEsf(i, j, k).y / (gradE_mag(i, j, k) + eps);
      const double nz = gradEsf(i, j, k).z / (gradE_mag(i, j, k) + eps);

      const double f = compute_eddington_factor(lambda(i, j, k), R(i, j, k));

      P_tensor(i, j, k).xx =
        (0.5 * (1 - f) + 0.5 * (3 * f - 1) * nx * nx) * Esf(i, j, k);
      P_tensor(i, j, k).xy = (0.5 * (3 * f - 1) * nx * ny) * Esf(i, j, k);
      P_tensor(i, j, k).xz = (0.5 * (3 * f - 1) * nx * nz) * Esf(i, j, k);

      P_tensor(i, j, k).yx = P_tensor(i, j, k).xy;
      P_tensor(i, j, k).yy =
        (0.5 * (1 - f) + 0.5 * (3 * f - 1) * ny * ny) * Esf(i, j, k);
      P_tensor(i, j, k).yz = (0.5 * (3 * f - 1) * ny * nz) * Esf(i, j, k);

      P_tensor(i, j, k).zx = P_tensor(i, j, k).xz;
      P_tensor(i, j, k).zy = P_tensor(i, j, k).yz;
      P_tensor(i, j, k).zz =
        (0.5 * (1 - f) + 0.5 * (3 * f - 1) * nz * nz) * Esf(i, j, k);
    };
  }
} // getTensorP

// Get the gradient of the radiaton energy density (E) using a 5-point stencil.
// The privilege for gradEsf_a should be `<wo, na>`, however it causes issues
// for the Legion tracing so I am using `<wo, ro>` for now.
template<std::size_t D>
void
getGradE(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, ro> Esf_a,
  typename field<vec<D>>::template accessor<wo, ro> gradEsf_a) {
  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto gradEsf = m.template mdcolex<is::cells>(gradEsf_a);

  if constexpr(D == 1) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());

    // Application of the 5-stencil central differencing:
    forall(i, (m.template cells<ax::x, dm::quantities>()), "getGradE") {
      gradEsf(i).x =
        (Esf(i - 2) - 8.0 * Esf(i - 1) + 8.0 * Esf(i + 1) - Esf(i + 2)) *
        one_over_12dx;
    }; // for
  }
  else if constexpr(D == 2) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());
    const double one_over_12dy = 1.0 / (12.0 * m.template delta<ax::y>());

    auto mdpolicy_qq = get_mdiota_policy(Esf,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "getGradE") {
      auto [j, i] = ji;

      // Application of the 5-stencil central differencing:
      gradEsf(i, j).x = (Esf(i - 2, j) - 8.0 * Esf(i - 1, j) +
                          8.0 * Esf(i + 1, j) - Esf(i + 2, j)) *
                        one_over_12dx;
      gradEsf(i, j).y = (Esf(i, j - 2) - 8.0 * Esf(i, j - 1) +
                          8.0 * Esf(i, j + 1) - Esf(i, j + 2)) *
                        one_over_12dy;
    }; // forall
  }
  else {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());
    const double one_over_12dy = 1.0 / (12.0 * m.template delta<ax::y>());
    const double one_over_12dz = 1.0 / (12.0 * m.template delta<ax::z>());

    auto mdpolicy_qqq = get_mdiota_policy(Esf,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "getGradE") {
      auto [k, j, i] = kji;

      // Application of the 5-stencil central differencing:
      gradEsf(i, j, k).x = (Esf(i - 2, j, k) - 8.0 * Esf(i - 1, j, k) +
                             8.0 * Esf(i + 1, j, k) - Esf(i + 2, j, k)) *
                           one_over_12dx;
      gradEsf(i, j, k).y = (Esf(i, j - 2, k) - 8.0 * Esf(i, j - 1, k) +
                             8.0 * Esf(i, j + 1, k) - Esf(i, j + 2, k)) *
                           one_over_12dy;
      gradEsf(i, j, k).z = (Esf(i, j, k - 2) - 8.0 * Esf(i, j, k - 1) +
                             8.0 * Esf(i, j, k + 1) - Esf(i, j, k + 2)) *
                           one_over_12dz;
    };
  }
} // getGradE

//
// Compute
//  1) magnitude of grad(E)
//  2) variable R (Eq 16 of Moens 2022 paper)
//  3) flux limiter function `lambda`
//
// Note : Lambda is only computed on the main grid (dm::quantities), then the
// outermost values are copied into ghost zones.
//
template<std::size_t D>
void
getLambda(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  field<double>::accessor<ro, na> Esf_a,
  typename field<vec<D>>::template accessor<ro, na> gradEsf_a,
  field<double>::accessor<rw, na> gradE_mag_a,
  field<double>::accessor<wo, na> R_a,
  field<double>::accessor<wo, na> lambda_a,
  single<double>::accessor<ro> kappa_a) {

  // TODO: applying boundary condition should be separated to new task

  auto r = m.template mdcolex<is::cells>(r_a);
  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto gradEsf = m.template mdcolex<is::cells>(gradEsf_a);
  auto gradE_mag = m.template mdcolex<is::cells>(gradE_mag_a);
  auto R = m.template mdcolex<is::cells>(R_a);
  auto lambda = m.template mdcolex<is::cells>(lambda_a);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "getLambda") {
      auto const kappa = *kappa_a;
      const double eps = 1.0e-30;

      gradE_mag(i) = std::abs(gradEsf(i).x);
      R(i) = gradE_mag(i) / (kappa * r(i) * Esf(i) + eps);
      lambda(i) = (2.0 + R(i)) / (6.0 + 3.0 * R(i) + R(i) * R(i));
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(Esf,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    forall(ji, mdpolicy_qq, "getLambda") {
      auto const kappa = *kappa_a;
      const double eps = 1.0e-30;

      auto [j, i] = ji;
      gradE_mag(i, j) = gradEsf(i, j).norm();
      R(i, j) = gradE_mag(i, j) / (kappa * r(i, j) * Esf(i, j) + eps);
      lambda(i, j) =
        (2.0 + R(i, j)) / (6.0 + 3.0 * R(i, j) + R(i, j) * R(i, j));
    };
  }
  else {
    const std::size_t xsize = m.template size<ax::x, dm::quantities>();
    const std::size_t ysize = m.template size<ax::y, dm::quantities>();
    const std::size_t zsize = m.template size<ax::z, dm::quantities>();
    auto mdpolicy_qqq = get_mdiota_policy(Esf,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    forall(kji, mdpolicy_qqq, "getLambda") {
      auto const kappa = *kappa_a;
      const double eps = 1.0e-30;

      auto [k, j, i] = kji;
      gradE_mag(i, j, k) = gradEsf(i, j, k).norm();
      R(i, j, k) =
        gradE_mag(i, j, k) / (kappa * r(i, j, k) * Esf(i, j, k) + eps);
      lambda(i, j, k) =
        (2.0 + R(i, j, k)) / (6.0 + 3.0 * R(i, j, k) + R(i, j, k) * R(i, j, k));
    };
  }
} // getLambda

// Get the radiation force using the FLD approximation
template<std::size_t D>
void
getRadForce(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> lambda_a,
  typename field<vec<D>>::template accessor<ro, na> gradEsf_a,
  typename field<vec<D>>::template accessor<wo, na> fr_a) {

  auto lambda = m.template mdcolex<is::cells>(lambda_a);
  auto gradEsf = m.template mdcolex<is::cells>(gradEsf_a);
  auto fr = m.template mdcolex<is::cells>(fr_a);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "getRadForce") {
      fr(i).x = -lambda(i) * gradEsf(i).x;
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(lambda,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "getRadForce") {
      auto [j, i] = ji;
      fr(i, j).x = -lambda(i, j) * gradEsf(i, j).x;
      fr(i, j).y = -lambda(i, j) * gradEsf(i, j).y;
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(lambda,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    forall(kji, mdpolicy_qqq, "getRadForce") {
      auto [k, j, i] = kji;
      fr(i, j, k).x = -lambda(i, j, k) * gradEsf(i, j, k).x;
      fr(i, j, k).y = -lambda(i, j, k) * gradEsf(i, j, k).y;
      fr(i, j, k).z = -lambda(i, j, k) * gradEsf(i, j, k).z;
    };
  }
} // getRadForce

// Radiation force and work terms. (Explicit source terms) (See Eq(21) in
// Moens2022)
template<std::size_t D>
void
explicitSourceUpdate(typename mesh<D>::template accessor<ro> m,
  // Primitive variables
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  // Terms required for computing radiation force and photon tiring
  typename field<vec<D>>::template accessor<ro, na> fr_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<ro,
    na> P_tensor_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<ro,
    na> gradV_a,
  // time derivative
  typename field<vec<D>>::template accessor<rw, na> dt_momentum_density_a,
  field<double>::accessor<rw, na> dt_total_energy_density_a,
  field<double>::accessor<rw, na> dt_radiation_energy_density_a) {

  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto fr = m.template mdcolex<is::cells>(fr_a);
  auto P_tensor = m.template mdcolex<is::cells>(P_tensor_a);
  auto gradV = m.template mdcolex<is::cells>(gradV_a);
  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);

  if constexpr(D == 1) {
    forall(
      i, (m.template cells<ax::x, dm::quantities>()), "explicitSourceUpdate") {

      // Explicitly updating conservative variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i) += fr(i);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i) += velocity(i).x * fr(i).x;

      // Subtracting the photon tiring term, (P::gradV), from the radiation
      // energy density in each cell. See Eq(34) in Moens2022.
      dt_radiation_energy_density(i) += -P_tensor(i).xx * gradV(i).xx;
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "explicitSourceUpdate") {
      auto [j, i] = ji;

      // Explicitly updating conservative variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j) += fr(i, j);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i, j) +=
        velocity(i, j).x * fr(i, j).x + velocity(i, j).y * fr(i, j).y;

      // Subtracting the photon tiring term, (P::gradV), from the radiation
      // energy density in each cell. See Eq(34) in Moens2022.
      dt_radiation_energy_density(i, j) +=
        -(P_tensor(i, j).xx * gradV(i, j).xx +
          P_tensor(i, j).xy * gradV(i, j).xy +
          P_tensor(i, j).yx * gradV(i, j).yx +
          P_tensor(i, j).yy * gradV(i, j).yy);
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "explicitSourceUpdate") {
      auto [k, j, i] = kji;

      // Explicitly updating conservative variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j, k) += fr(i, j, k);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x * fr(i).x
      dt_total_energy_density(i, j, k) += velocity(i, j, k).x * fr(i, j, k).x +
                                          velocity(i, j, k).y * fr(i, j, k).y +
                                          velocity(i, j, k).z * fr(i, j, k).z;

      // Subtracting the photon tiring term, (P::gradV), from the radiation
      // energy density in each cell. See Eq(34) in Moens et al. 2022.
      dt_radiation_energy_density(i, j, k) +=
        -(P_tensor(i, j, k).xx * gradV(i, j, k).xx +
          P_tensor(i, j, k).xy * gradV(i, j, k).xy +
          P_tensor(i, j, k).xz * gradV(i, j, k).xz +
          P_tensor(i, j, k).yx * gradV(i, j, k).yx +
          P_tensor(i, j, k).yy * gradV(i, j, k).yy +
          P_tensor(i, j, k).yz * gradV(i, j, k).yz +
          P_tensor(i, j, k).zx * gradV(i, j, k).zx +
          P_tensor(i, j, k).zy * gradV(i, j, k).zy +
          P_tensor(i, j, k).zz * gradV(i, j, k).zz);
    }; // forall
  }
} // explicitSourceUpdate

// Compute the diffusion coefficients D
template<std::size_t D>
void
getDiff(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  field<double>::accessor<ro, na> lambda_a,
  field<double>::accessor<wo, na> Diff_a,
  single<double>::accessor<ro> kappa_a) {

  auto r = m.template mdcolex<is::cells>(r_a);
  auto lambda = m.template mdcolex<is::cells>(lambda_a);
  auto Diff = m.template mdcolex<is::cells>(Diff_a);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "getDiff") {
      auto const kappa = *kappa_a;
      const double clight = hard::constants::cgs::speed_of_light;
      Diff(i) = clight * lambda(i) / (kappa * r(i));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_aa = get_mdiota_policy(r,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_aa, "getDiff") {
      auto const kappa = *kappa_a;
      const double clight = hard::constants::cgs::speed_of_light;
      auto [j, i] = ji;
      Diff(i, j) = clight * lambda(i, j) / (kappa * r(i, j));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_aaa = get_mdiota_policy(r,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_aaa, "getDiff") {
      auto const kappa = *kappa_a;
      const double clight = hard::constants::cgs::speed_of_light;
      auto [k, j, i] = kji;
      Diff(i, j, k) = clight * lambda(i, j, k) / (kappa * r(i, j, k));
    };
  } // if
} // getDiff

// Compute the radiation diffusion coefficients D on faces
template<std::size_t D>
void
diffusion_init(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> Diff_a,
  typename field<double>::template accessor<wo, na> Df_xa,
  typename field<double>::template accessor<wo, na> Df_ya,
  typename field<double>::template accessor<wo, na> Df_za) {

  auto Diff = m.template mdcolex<is::cells>(Diff_a);
  auto Df_x = m.template mdcolex<is::cells>(Df_xa);
  auto Df_y = m.template mdcolex<is::cells>(Df_ya);
  auto Df_z = m.template mdcolex<is::cells>(Df_za);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::corrector>()), "diffusion_init") {
      Df_x(i) = 2 * Diff(i) * Diff(i - 1) / (Diff(i) + Diff(i - 1));
    };
  }
  else if constexpr(D == 2) {

    auto mdpolicy_cc = get_mdiota_policy(Diff,
      m.template cells<ax::y, dm::corrector>(),
      m.template cells<ax::x, dm::corrector>());

    forall(ji, mdpolicy_cc, "explicitSourceUpdate") {
      auto [j, i] = ji;
      Df_x(i, j) =
        2 * Diff(i, j) * Diff(i - 1, j) / (Diff(i, j) + Diff(i - 1, j));
      Df_y(i, j) =
        2 * Diff(i, j) * Diff(i, j - 1) / (Diff(i, j) + Diff(i, j - 1));
    }; // for
  }
  else /* D == 3 */ {
    auto mdpolicy_ccc = get_mdiota_policy(Diff,
      m.template cells<ax::z, dm::corrector>(),
      m.template cells<ax::y, dm::corrector>(),
      m.template cells<ax::x, dm::corrector>());

    forall(kji, mdpolicy_ccc, "diffusion_init") {
      auto [k, j, i] = kji;
      Df_x(i, j, k) = 2 * Diff(i, j, k) * Diff(i - 1, j, k) /
                      (Diff(i, j, k) + Diff(i - 1, j, k));
      Df_y(i, j, k) = 2 * Diff(i, j, k) * Diff(i, j - 1, k) /
                      (Diff(i, j, k) + Diff(i, j - 1, k));
      Df_z(i, j, k) = 2 * Diff(i, j, k) * Diff(i, j, k - 1) /
                      (Diff(i, j, k) + Diff(i, j, k - 1));
    };
  } // if
} // diffusion_init

template<std::size_t D>
void
const_init(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<wo, na> f_a,
  double w) {

  auto f = m.template mdcolex<is::cells>(f_a);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "const_init") {
      f(i) = w;
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_aa = get_mdiota_policy(f,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_aa, "const_init") {
      auto [j, i] = ji;
      f(i, j) = w;
    };
  }
  else /* D == 3 */ {
    // TODO: Loop collapse for these tiny tasks not a good idea
    auto mdpolicy_aaa = get_mdiota_policy(f,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_aaa, "const_init") {
      auto [k, j, i] = kji;
      f(i, j, k) = w;
    };
  } // if
} // const_init

template<std::size_t D>
void
initialize_Ef(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, na> Erad_a,
  typename field<double>::template accessor<wo, na> Ef_a) {

  auto Erad = m.template mdcolex<is::cells>(Erad_a);
  auto Ef = m.template mdcolex<is::cells>(Ef_a);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "initialize_Ef") {
      Ef(i) = Erad(i);
    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(Ef,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "initialize_Ef") {
      auto [j, i] = ji;
      Ef(i, j) = Erad(i, j);
    }; // forall
  }
  else if constexpr(D == 3) {
    auto mdpolicy_qqq = get_mdiota_policy(Ef,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "initialize_Ef") {
      auto [k, j, i] = kji;
      Ef(i, j, k) = Erad(i, j, k);
    }; // forall
  }
} // initialize_Ef

template<std::size_t D>
void
store_du_dt_implicit_from_diffusion(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, na> Esf_a,
  typename field<double>::template accessor<ro, na> Ef_a,
  single<double>::accessor<ro> dtw_a,
  typename field<double>::template accessor<rw, na>
    dt_radiation_energy_density_implicit_a) {

  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto Ef = m.template mdcolex<is::cells>(Ef_a);
  auto dt_radiation_energy_density_implicit =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_implicit_a);

  if constexpr(D == 1) {
    forall(i,
      (m.template cells<ax::x, dm::quantities>()),
      "store_du_dt_implicit_from_diffusion") {
      auto & dt_weighted = *dtw_a;
      const double one_over_aii_dt{1.0 / dt_weighted};
      dt_radiation_energy_density_implicit(i) +=
        (Esf(i) - Ef(i)) * one_over_aii_dt;
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(Ef,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "initialize_Ef") {
      auto [j, i] = ji;
      auto & dt_weighted = *dtw_a;
      const double one_over_aii_dt{1.0 / dt_weighted};
      dt_radiation_energy_density_implicit(i, j) +=
        (Esf(i, j) - Ef(i, j)) * one_over_aii_dt;
    }; // for
  }
  else if constexpr(D == 3) {
    auto mdpolicy_qqq = get_mdiota_policy(Ef,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "initialize_Ef") {
      auto [k, j, i] = kji;
      auto & dt_weighted = *dtw_a;
      const double one_over_aii_dt{1.0 / dt_weighted};
      dt_radiation_energy_density_implicit(i, j, k) +=
        (Esf(i, j, k) - Ef(i, j, k)) * one_over_aii_dt;
    }; // for
  }
} // copy_esf_into_ef

template<std::size_t D>
void
stencil_init(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> Df_xa,
  typename field<double>::template accessor<ro, ro> Df_ya,
  typename field<double>::template accessor<ro, ro> Df_za,
  typename field<stencil<D>>::template accessor<wo, na> Ew_a,
  single<double>::accessor<ro> dt_a) {
  // TODO: Stencil, Ew ghosts can be `na` (?)

  auto Df_x = m.template mdcolex<is::cells>(Df_xa);
  auto Df_y = m.template mdcolex<is::cells>(Df_ya);
  auto Df_z = m.template mdcolex<is::cells>(Df_za);
  auto Ew = m.template mdcolex<is::cells>(Ew_a);

  if constexpr(D == 1) {
    const double dx{m.template delta<ax::x>()};

    forall(i, (m.template cells<ax::x, dm::corrector>()), "stencil_init") {
      auto const dt = *dt_a;
      const double wx{dt / pow(dx, 2)};
      Ew(i)[dirs::c] = 1.0 + wx * (Df_x(i + 1) + Df_x(i));
      Ew(i)[dirs::w] = wx * Df_x(i);
    };
  }
  else if constexpr(D == 2) {
    const double dx{m.template delta<ax::x>()};
    const double dy{m.template delta<ax::y>()};

    auto mdpolicy_qq = get_mdiota_policy(Df_x,
      m.template cells<ax::y, dm::corrector>(),
      m.template cells<ax::x, dm::corrector>());

    forall(ji, mdpolicy_qq, "explicitSourceUpdate") {
      auto [j, i] = ji;
      auto const dt = *dt_a;
      const double wx{dt / pow(dx, 2)};
      const double wy{dt / pow(dy, 2)};
      Ew(i, j)[dirs::c] = 1.0 + (wx * (Df_x(i + 1, j) + Df_x(i, j)) +
                                  wy * (Df_y(i, j + 1) + Df_y(i, j)));
      Ew(i, j)[dirs::w] = wx * Df_x(i, j);
      Ew(i, j)[dirs::s] = wy * Df_y(i, j);
    }; // for
  }
  else /* D == 3 */ {
    const double dx{m.template delta<ax::x>()};
    const double dy{m.template delta<ax::y>()};
    const double dz{m.template delta<ax::z>()};

    auto mdpolicy_qqq = get_mdiota_policy(Df_x,
      m.template cells<ax::z, dm::corrector>(),
      m.template cells<ax::y, dm::corrector>(),
      m.template cells<ax::x, dm::corrector>());

    forall(kji, mdpolicy_qqq, "stencil_init") {
      auto [k, j, i] = kji;
      auto const dt = *dt_a;
      const double wx{dt / pow(dx, 2)};
      const double wy{dt / pow(dy, 2)};
      const double wz{dt / pow(dz, 2)};
      Ew(i, j, k)[dirs::c] = 1.0 + (wx * (Df_x(i + 1, j, k) + Df_x(i, j, k)) +
                                     wy * (Df_y(i, j + 1, k) + Df_y(i, j, k)) +
                                     wz * (Df_z(i, j, k + 1) + Df_z(i, j, k)));
      Ew(i, j, k)[dirs::w] = wx * Df_x(i, j, k);
      Ew(i, j, k)[dirs::s] = wy * Df_y(i, j, k);
      Ew(i, j, k)[dirs::d] = wz * Df_z(i, j, k);
    };
  } // if
} // stencil_init

template<std::size_t D>
void
apply_BC(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, wo> E_a) {

  auto E = m.template mdcolex<is::cells>(E_a);
  // TODO: Other D?
  if constexpr(D == 1) {
    const std::size_t i = m.template size<ax::x, dm::all>();
    // periodic boundary for 1D
    E(0) = E(i - 4);
    E(1) = E(i - 3);
    E(i - 2) = E(2);
    E(i - 1) = E(3);
  }
} // apply_BC

template<std::size_t D>
void
full_weighting(typename mesh<D>::template accessor<ro> mf,
  typename mesh<D>::template accessor<ro> mc,
  typename field<double>::template accessor<ro, ro> rfa,
  typename field<double>::template accessor<wo, ro> fca) {
  // TODO: fca could be <wo, na>, since only writing quantities

  auto rf = mf.template mdcolex<is::cells>(rfa);
  auto fc = mc.template mdcolex<is::cells>(fca);

  if constexpr(D == 1) {
    forall(i, (mc.template cells<ax::x, dm::quantities>()), "full_weighting") {
      auto fi = 2 * i;
      fc(i) = 0.125 * (3.0 * (rf(fi - 2) + rf(fi - 1)) + rf(fi - 3) + rf(fi));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(fc,
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "full_weighting") {
      auto [j, i] = ji;
      auto fj = 2 * j;
      auto fi = 2 * i;
      fc(i, j) =
        0.015625 *
        (9.0 * (rf(fi - 2, fj - 2) + rf(fi - 1, fj - 2) + rf(fi - 2, fj - 1) +
                 rf(fi - 1, fj - 1)) +

          3.0 * (rf(fi - 3, fj - 2) + rf(fi - 3, fj - 1) + rf(fi, fj - 2) +
                  rf(fi, fj - 1) + rf(fi - 2, fj - 3) + rf(fi - 1, fj - 3) +
                  rf(fi - 2, fj) + rf(fi - 1, fj)) +

          rf(fi, fj - 3) + rf(fi - 3, fj - 3) + rf(fi - 3, fj) + rf(fi, fj));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(fc,
      mc.template cells<ax::z, dm::quantities>(),
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "full_weighting") {
      auto [k, j, i] = kji;
      auto fk = 2 * k;
      auto fj = 2 * j;
      auto fi = 2 * i;

      fc(i, j, k) =
        0.001953125 *
        (27.0 * (rf(fi - 2, fj - 2, fk - 2) + rf(fi - 2, fj - 2, fk - 1) +
                  rf(fi - 2, fj - 1, fk - 2) + rf(fi - 1, fj - 2, fk - 2) +
                  rf(fi - 2, fj - 1, fk - 1) + rf(fi - 1, fj - 2, fk - 1) +
                  rf(fi - 1, fj - 1, fk - 2) + rf(fi - 1, fj - 1, fk - 1)) +

          9.0 * (rf(fi - 2, fj - 2, fk - 3) + rf(fi - 2, fj - 3, fk - 2) +
                  rf(fi - 3, fj - 2, fk - 2) + rf(fi - 2, fj - 1, fk - 3) +
                  rf(fi - 2, fj - 3, fk - 1) + rf(fi - 3, fj - 2, fk - 1) +
                  rf(fi - 1, fj - 2, fk - 3) + rf(fi - 1, fj - 3, fk - 2) +
                  rf(fi - 3, fj - 1, fk - 2) + rf(fi - 1, fj - 1, fk - 3) +
                  rf(fi - 1, fj - 3, fk - 1) + rf(fi - 3, fj - 1, fk - 1) +

                  rf(fi - 2, fj - 2, fk) + rf(fi - 2, fj, fk - 2) +
                  rf(fi, fj - 2, fk - 2) + rf(fi - 2, fj - 1, fk) +
                  rf(fi - 2, fj, fk - 1) + rf(fi, fj - 2, fk - 1) +
                  rf(fi - 1, fj - 2, fk) + rf(fi - 1, fj, fk - 2) +
                  rf(fi, fj - 1, fk - 2) + rf(fi - 1, fj - 1, fk) +
                  rf(fi - 1, fj, fk - 1) + rf(fi, fj - 1, fk - 1)) +

          3.0 *
            (rf(fi - 2, fj - 3, fk - 3) + rf(fi - 3, fj - 2, fk - 3) +
              rf(fi - 3, fj - 3, fk - 2) + rf(fi - 1, fj - 3, fk - 3) +
              rf(fi - 3, fj - 1, fk - 3) + rf(fi - 3, fj - 3, fk - 1) +

              rf(fi - 2, fj, fk - 3) + rf(fi, fj - 2, fk - 3) +
              rf(fi, fj - 3, fk - 2) + rf(fi - 1, fj, fk - 3) +
              rf(fi, fj - 1, fk - 3) + rf(fi, fj - 3, fk - 1) +

              rf(fi - 2, fj - 3, fk) + rf(fi - 3, fj - 2, fk) +
              rf(fi - 3, fj, fk - 2) + rf(fi - 1, fj - 3, fk) +
              rf(fi - 3, fj - 1, fk) + rf(fi - 3, fj, fk - 1) +

              rf(fi - 2, fj, fk) + rf(fi, fj - 2, fk) + rf(fi, fj, fk - 2) +
              rf(fi - 1, fj, fk) + rf(fi, fj - 1, fk) + rf(fi, fj, fk - 1)) +

          rf(fi - 3, fj - 3, fk - 3) + rf(fi - 3, fj - 3, fk) +
          rf(fi - 3, fj, fk - 3) + rf(fi, fj - 3, fk - 3) + rf(fi - 3, fj, fk) +
          rf(fi, fj - 3, fk) + rf(fi, fj, fk - 3) + rf(fi, fj, fk));
    }; // for
  } // if
} // full_weighting

template<std::size_t D>
void
nlinear_interpolation(typename mesh<D>::template accessor<ro> mc,
  typename mesh<D>::template accessor<ro> mf,
  typename field<double>::template accessor<ro, ro> cfa,
  typename field<double>::template accessor<wo, ro> ffa) {
  // TODO: As above, ffa could be <wo, na>, since only writing quantities

  auto cf = mc.template mdcolex<is::cells>(cfa);
  auto ff = mf.template mdcolex<is::cells>(ffa);

  if constexpr(D == 1) {
    forall(i,
      (mf.template cells<ax::x, dm::quantities>()),
      "nlinear_interpolation") {
      auto ci = i / 2 + 1;
      auto di = +1;
      if((i + 1) % 2) {
        auto di = -1;
      }
      ff(i) = 0.25 * (3 * cf(ci) + cf(ci + di));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(cf,
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "full_weighting") {
      auto [j, i] = ji;
      auto cj = j / 2 + 1;
      auto dj = 1;
      if((j + 1) % 2) {
        auto dj = -1;
      }
      auto ci = i / 2 + 1;
      auto di = 1;
      if((i + 1) % 2) {
        auto di = -1;
      }
      ff(i, j) = 0.0625 * (9 * cf(ci, cj) + 3 * cf(ci + di, cj) +
                            3 * cf(ci, cj + dj) + cf(ci + di, cj + dj));
    }; // for
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(cf,
      mc.template cells<ax::z, dm::quantities>(),
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "full_weighting") {
      auto [k, j, i] = kji;
      auto ck = k / 2 + 1;
      auto dk = 1;
      if((k + 1) % 2) {
        auto dk = -1;
      }
      auto cj = j / 2 + 1;
      auto dj = 1;
      if((j + 1) % 2) {
        auto dj = -1;
      }
      auto ci = i / 2 + 1;
      auto di = 1;
      if((i + 1) % 2) {
        auto di = -1;
      }
      ff(i, j, k) =
        0.015625 *
        (27 * cf(ci, cj, ck) + 9 * cf(ci + di, cj, ck) +
          9 * cf(ci, cj + dj, ck) + 9 * cf(ci, cj, ck + dk) +
          3 * cf(ci + di, cj + dj, ck) + 3 * cf(ci + di, cj, ck + dk) +
          3 * cf(ci, cj + dj, ck + dk) + cf(ci + di, cj + dj, ck + dk));
    };
  } // if
} // nlinear_interpolation

template<std::size_t D>
void
damped_jacobi(typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<rw, ro> ua_new,
  typename field<double>::template accessor<ro, ro> ua_old,
  typename field<double>::template accessor<ro, ro> fa,
  double omega) {

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u_new = m.template mdcolex<is::cells>(ua_new);
  auto u_old = m.template mdcolex<is::cells>(ua_old);
  auto f = m.template mdcolex<is::cells>(fa);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "damped_jacobi") {
      const double z = (Ew(i)[dirs::w] * u_old(i - 1) +
                         Ew(i + 1)[dirs::w] * u_old(i + 1) + f(i)) /
                       Ew(i)[dirs::c];

      u_new(i) = u_old(i) + omega * (z - u_old(i));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u_new,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "damped_jacobi") {
      auto [j, i] = ji;
      const double z = (Ew(i, j)[dirs::w] * u_old(i - 1, j) +
                         Ew(i + 1, j)[dirs::w] * u_old(i + 1, j) +
                         Ew(i, j)[dirs::s] * u_old(i, j - 1) +
                         Ew(i, j + 1)[dirs::s] * u_old(i, j + 1) + f(i, j)) /
                       Ew(i, j)[dirs::c];

      u_new(i, j) = u_old(i, j) + omega * (z - u_old(i, j));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u_new,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "damped_jacobi") {
      auto [k, j, i] = kji;
      const double z =
        (Ew(i, j, k)[dirs::w] * u_old(i - 1, j, k) +
          Ew(i + 1, j, k)[dirs::w] * u_old(i + 1, j, k) +
          Ew(i, j, k)[dirs::s] * u_old(i, j - 1, k) +
          Ew(i, j + 1, k)[dirs::s] * u_old(i, j + 1, k) +
          Ew(i, j, k)[dirs::d] * u_old(i, j, k - 1) +
          Ew(i, j, k + 1)[dirs::d] * u_old(i, j, k + 1) + f(i, j, k)) /
        Ew(i, j, k)[dirs::c];

      u_new(i, j, k) = u_old(i, j, k) + omega * (z - u_old(i, j, k));
    }; // forall
  } // if
} // damped_jacobi

// Gauss-Seidel
template<std::size_t D>
void
gauss_seidel(typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<rw, ro> ua,
  typename field<double>::template accessor<ro, ro> fa) {

  // TODO: Parallelize Gauss seidel

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u = m.template mdcolex<is::cells>(ua);
  auto f = m.template mdcolex<is::cells>(fa);

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      u(i) =
        (Ew(i)[dirs::w] * u(i - 1) + Ew(i + 1)[dirs::w] * u(i + 1) + f(i)) /
        Ew(i)[dirs::c];

    } // for
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        u(i, j) = (Ew(i, j)[dirs::w] * u(i - 1, j) +
                    Ew(i + 1, j)[dirs::w] * u(i + 1, j) +
                    Ew(i, j)[dirs::s] * u(i, j - 1) +
                    Ew(i, j + 1)[dirs::s] * u(i, j + 1) + f(i, j)) /
                  Ew(i, j)[dirs::c];

      } // for
    } // for
  }
  else /* D == 3 */ {
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          u(i, j, k) =
            (Ew(i, j, k)[dirs::w] * u(i - 1, j, k) +
              Ew(i + 1, j, k)[dirs::w] * u(i + 1, j, k) +
              Ew(i, j, k)[dirs::s] * u(i, j - 1, k) +
              Ew(i, j + 1, k)[dirs::s] * u(i, j + 1, k) +
              Ew(i, j, k)[dirs::d] * u(i, j, k - 1) +
              Ew(i, j, k + 1)[dirs::d] * u(i, j, k + 1) + f(i, j, k)) /
            Ew(i, j, k)[dirs::c];

        } // for
      } // for
    } // for
  } // if
} // gauss_seidel

// Red Gauss-Seidel
template<std::size_t D>
void
red_gauss_seidel(typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<rw, ro> ua,
  typename field<double>::template accessor<ro, ro> fa) {

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u = m.template mdcolex<is::cells>(ua);
  auto f = m.template mdcolex<is::cells>(fa);

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {

      unsigned int id_i = i;

      if((id_i % 2) == 0) {
        u(i) =
          (Ew(i)[dirs::w] * u(i - 1) + Ew(i + 1)[dirs::w] * u(i + 1) + f(i)) /
          Ew(i)[dirs::c];
      }

    } // for
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {

        unsigned int id_i = i;
        unsigned int id_j = j;

        if(((id_i + id_j) % 2) == 0) {
          u(i, j) = (Ew(i, j)[dirs::w] * u(i - 1, j) +
                      Ew(i + 1, j)[dirs::w] * u(i + 1, j) +
                      Ew(i, j)[dirs::s] * u(i, j - 1) +
                      Ew(i, j + 1)[dirs::s] * u(i, j + 1) + f(i, j)) /
                    Ew(i, j)[dirs::c];
        }

      } // for
    } // for
  }
  else /* D == 3 */ {
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {

          unsigned int id_i = i;
          unsigned int id_j = j;
          unsigned int id_k = k;

          if(((id_i + id_j + id_k) % 2) == 0) {
            u(i, j, k) =
              (Ew(i, j, k)[dirs::w] * u(i - 1, j, k) +
                Ew(i + 1, j, k)[dirs::w] * u(i + 1, j, k) +
                Ew(i, j, k)[dirs::s] * u(i, j - 1, k) +
                Ew(i, j + 1, k)[dirs::s] * u(i, j + 1, k) +
                Ew(i, j, k)[dirs::d] * u(i, j, k - 1) +
                Ew(i, j, k + 1)[dirs::d] * u(i, j, k + 1) + f(i, j, k)) /
              Ew(i, j, k)[dirs::c];
          }

        } // for
      } // for
    } // for
  } // if
} // red_gauss_seidel

// Black Gauss-Seidel
template<std::size_t D>
void
black_gauss_seidel(typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<rw, ro> ua,
  typename field<double>::template accessor<ro, ro> fa) {

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u = m.template mdcolex<is::cells>(ua);
  auto f = m.template mdcolex<is::cells>(fa);

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {

      unsigned int id_i = i;

      if((id_i % 2) == 1) {
        u(i) =
          (Ew(i)[dirs::w] * u(i - 1) + Ew(i + 1)[dirs::w] * u(i + 1) + f(i)) /
          Ew(i)[dirs::c];
      }

    } // for
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {

        unsigned int id_i = i;
        unsigned int id_j = j;

        if(((id_i + id_j) % 2) == 1) {
          u(i, j) = (Ew(i, j)[dirs::w] * u(i - 1, j) +
                      Ew(i + 1, j)[dirs::w] * u(i + 1, j) +
                      Ew(i, j)[dirs::s] * u(i, j - 1) +
                      Ew(i, j + 1)[dirs::s] * u(i, j + 1) + f(i, j)) /
                    Ew(i, j)[dirs::c];
        }
      } // for
    } // for
  }
  else /* D == 3 */ {
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {

          unsigned int id_i = i;
          unsigned int id_j = j;
          unsigned int id_k = k;

          if(((id_i + id_j + id_k) % 2) == 1) {
            u(i, j, k) =
              (Ew(i, j, k)[dirs::w] * u(i - 1, j, k) +
                Ew(i + 1, j, k)[dirs::w] * u(i + 1, j, k) +
                Ew(i, j, k)[dirs::s] * u(i, j - 1, k) +
                Ew(i, j + 1, k)[dirs::s] * u(i, j + 1, k) +
                Ew(i, j, k)[dirs::d] * u(i, j, k - 1) +
                Ew(i, j, k + 1)[dirs::d] * u(i, j, k + 1) + f(i, j, k)) /
              Ew(i, j, k)[dirs::c];
          }
        } // for
      } // for
    } // for
  } // if
} // black_gauss_seidel

template<std::size_t D>
void
residual(typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<ro, ro> ua,
  typename field<double>::template accessor<ro, ro> fa,
  typename field<double>::template accessor<wo, ro> ra) {
  // TODO: Not using the ghost cells of fa here, are we?

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u = m.template mdcolex<is::cells>(ua);
  auto f = m.template mdcolex<is::cells>(fa);
  auto r = m.template mdcolex<is::cells>(ra);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "residual") {
      r(i) = f(i) - (Ew(i)[dirs::c] * u(i) - Ew(i)[dirs::w] * u(i - 1) -
                      Ew(i + 1)[dirs::w] * u(i + 1));
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "residual") {
      auto [j, i] = ji;
      r(i, j) = f(i, j) -
                (Ew(i, j)[dirs::c] * u(i, j) - Ew(i, j)[dirs::w] * u(i - 1, j) -
                  Ew(i + 1, j)[dirs::w] * u(i + 1, j) -
                  Ew(i, j)[dirs::s] * u(i, j - 1) -
                  Ew(i, j + 1)[dirs::s] * u(i, j + 1));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "residual") {
      auto [k, j, i] = kji;
      r(i, j, k) = f(i, j, k) - (Ew(i, j, k)[dirs::c] * u(i, j, k) -
                                  Ew(i, j, k)[dirs::w] * u(i - 1, j, k) -
                                  Ew(i + 1, j, k)[dirs::w] * u(i + 1, j, k) -
                                  Ew(i, j, k)[dirs::s] * u(i, j - 1, k) -
                                  Ew(i, j + 1, k)[dirs::s] * u(i, j + 1, k) -
                                  Ew(i, j, k)[dirs::d] * u(i, j, k - 1) -
                                  Ew(i, j, k + 1)[dirs::d] * u(i, j, k + 1));
    }; // forall
  } // if
} // residual

template<std::size_t D>
void
correction(typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<rw, ro> ua,
  typename field<double>::template accessor<wo, ro> ea) {
  // TODO: Looks like ea should be <ro, na>

  auto u = m.template mdcolex<is::cells>(ua);
  auto e = m.template mdcolex<is::cells>(ea);

  if constexpr(D == 1) {
    forall(i, (m.template cells<ax::x, dm::quantities>()), "correction") {
      u(i) += e(i);
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(ji, mdpolicy_qq, "full_weighting") {
      auto [j, i] = ji;
      u(i, j) += e(i, j);
    }; // for
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    forall(kji, mdpolicy_qqq, "full_weighting") {
      auto [k, j, i] = kji;
      u(i, j, k) += e(i, j, k);
    }; // for
  } // if
} // correction

} // namespace hard::task::rad

#endif // HARD_TASKS_RAD_HH
