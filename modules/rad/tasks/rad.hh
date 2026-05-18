#ifndef HARD_MODULE_RAD_TASKS_RAD_HH
#define HARD_MODULE_RAD_TASKS_RAD_HH

#include "../constants.hh"
#include "utils.hh"
#include <../modules/spec/utils.hh>

namespace hard::tasks::rad {

// Castro/Source/radiation/fluxlimiter.H/Edd_factor
FLECSI_INLINE_TARGET double
AFLDEddFactor(double lam,
  std::size_t limiter_id,
  std::size_t closure_id) noexcept {
  double f = 0.0;

  switch(closure_id) {
    case 0: // f = λ
      f = lam;
      break;
    case 1: // f = 1/3
      f = 1.0 / 3.0;
      break;
    case 2: // f = 1 - 2λ
      f = 1.0 - 2.0 * lam;
      break;
    case 3: // f = λ + (λR)^2 with different limiters
      // lambda + (lambda*R)**2
      switch(limiter_id) {
        case 0:
          f = 1.0 / 3.0;
          break;
        case 1: {
          double temp =
            0.5 * std::max(0.0, 1.0 - 3.0 * lam) +
            std::sqrt(std::max(0.0, (1.0 - 3.0 * lam) * (1.0 + 5.0 * lam)));
          f = lam + temp * temp; // approximate LP, [123]
          break;
        }
        case 2:
          f = 1.0 - 5.0 * lam + 9.0 * lam * lam; // Bruenn, 1[123]
          break;
        case 3:
          f = 1.0 + lam - 9.0 * lam * lam; // Larsen's square root, 2[123]
          break;
        case 4:
          if(lam > 2.0 / 9.0) // Minerbo
            f = 1.0 / 3.0;
          else
            f = 1.0 + 3.0 * lam - 2.0 * std::sqrt(2.0 * lam);
          break;
        default:
          assert(false && "Invalid Limiter ID (Closure 3)");
          return -1.0;
      }
      break;
    case 4: // f = 1/3 + (2/3)*(λR)^2 with different limiters
      switch(limiter_id) {
        case 0:
          f = 1.0 / 3.0; // no limiter
          break;
        case 1: {
          double temp =
            std::max(0.0, 1.0 - 3.0 * lam) +
            std::sqrt(std::max(0.0, (1.0 - 3.0 * lam) * (1.0 + 5.0 * lam)));
          f = 1.0 / 3.0 + (temp * temp / 6.0); // approximate LP, [123]
          break;
        }
        case 2:
          f = 1.0 / 3.0 +
              2.0 * (1.0 - 6.0 * lam + 9.0 * lam * lam) / 3.0; // Bruenn, 1[123]
          break;
        case 3:
          f = 1.0 / 3.0 + 2.0 * (1.0 - 9.0 * lam * lam) /
                            3.0; // Larsen's square root, 2[123]
          break;
        case 4:
          if(lam > 2.0 / 9.0)
            f = 5.0 / 9.0 - (2.0 * lam / 3.0); // Minerbo
          else
            f = 1.0 / 3.0 +
                (2.0 * (1.0 + 2.0 * lam - 2.0 * std::sqrt(2.0 * lam)) / 3.0);
          break;
        default:
          assert(false && "Invalid Limiter ID (Closure 4)");
          return -1.0;
      }
      break;
    default:
      assert(false && "Invalid Closure ID");
      return -1.0;
  }

  return f;
}

// Castro/Source/radiation/rad_util.H/FLDlambda
FLECSI_INLINE_TARGET double
AFLDlambda(double r, std::size_t limiter_id) noexcept {
  double l = 0.0;

  switch(limiter_id) {
    case 0:
      l = 1.0 / 3.0; // no limiter
      break;
    case 1:
      l = (2.0 + r) / (6.0 + 3.0 * r + r * r); // approximate LP
      break;
    case 2:
      l = 1.0 / (3.0 + r); // Bruenn
      break;
    case 3:
      l = 1.0 / std::sqrt(9.0 + r * r); // Larsen's square root
      break;
    case 4:
      if(r < 1.5)
        l = 2.0 / (3.0 + std::sqrt(9.0 + 12.0 * r * r)); // Minerbo
      else
        l = 1.0 / (1.0 + r + std::sqrt(1.0 + 2.0 * r));
      break;
    default:
      assert(false && "Invalid Limiter ID");
      return -1.0;
  }

  return l;
}

using hard::tasks::util::get_mdiota_policy;

// Get the gradient of the velocity using a 5-point stencil.
template<std::size_t D>
void
getGradV(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    na> gradV_a,
  typename field<vec<D>>::template accessor<ro, ro> u_a) noexcept {

  auto u = m.template mdcolex<is::cells>(u_a);
  auto grad_v = m.template mdcolex<is::cells>(gradV_a);

  if constexpr(D == 1) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());

    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      grad_v(i).xx = (u(i - 2).x() - 8.0 * u(i - 1).x() + 8.0 * u(i + 1).x() -
                       u(i + 2).x()) *
                     one_over_12dx;
    };
  }
  else if constexpr(D == 2) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());
    const double one_over_12dy = 1.0 / (12.0 * m.template delta<ax::y>());

    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      grad_v(i, j).xx = (u(i - 2, j).x() - 8.0 * u(i - 1, j).x() +
                          8.0 * u(i + 1, j).x() - u(i + 2, j).x()) *
                        one_over_12dx;
      grad_v(i, j).xy = (u(i, j - 2).x() - 8.0 * u(i, j - 1).x() +
                          8.0 * u(i, j + 1).x() - u(i, j + 2).x()) *
                        one_over_12dy;

      grad_v(i, j).yx = (u(i - 2, j).y() - 8.0 * u(i - 1, j).y() +
                          8.0 * u(i + 1, j).y() - u(i + 2, j).y()) *
                        one_over_12dx;
      grad_v(i, j).yy = (u(i, j - 2).y() - 8.0 * u(i, j - 1).y() +
                          8.0 * u(i, j + 1).y() - u(i, j + 2).y()) *
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

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      grad_v(i, j, k).xx = (u(i - 2, j, k).x() - 8.0 * u(i - 1, j, k).x() +
                             8.0 * u(i + 1, j, k).x() - u(i + 2, j, k).x()) *
                           one_over_12dx;
      grad_v(i, j, k).xy = (u(i, j - 2, k).x() - 8.0 * u(i, j - 1, k).x() +
                             8.0 * u(i, j + 1, k).x() - u(i, j + 2, k).x()) *
                           one_over_12dy;
      grad_v(i, j, k).xz = (u(i, j, k - 2).x() - 8.0 * u(i, j, k - 1).x() +
                             8.0 * u(i, j, k + 1).x() - u(i, j, k + 2).x()) *
                           one_over_12dz;

      grad_v(i, j, k).yx = (u(i - 2, j, k).y() - 8.0 * u(i - 1, j, k).y() +
                             8.0 * u(i + 1, j, k).y() - u(i + 2, j, k).y()) *
                           one_over_12dx;
      grad_v(i, j, k).yy = (u(i, j - 2, k).y() - 8.0 * u(i, j - 1, k).y() +
                             8.0 * u(i, j + 1, k).y() - u(i, j + 2, k).y()) *
                           one_over_12dy;
      grad_v(i, j, k).yz = (u(i, j, k - 2).y() - 8.0 * u(i, j, k - 1).y() +
                             8.0 * u(i, j, k + 1).y() - u(i, j, k + 2).y()) *
                           one_over_12dz;

      grad_v(i, j, k).zx = (u(i - 2, j, k).z() - 8.0 * u(i - 1, j, k).z() +
                             8.0 * u(i + 1, j, k).z() - u(i + 2, j, k).z()) *
                           one_over_12dx;
      grad_v(i, j, k).zy = (u(i, j - 2, k).z() - 8.0 * u(i, j - 1, k).z() +
                             8.0 * u(i, j + 1, k).z() - u(i, j + 2, k).z()) *
                           one_over_12dy;
      grad_v(i, j, k).zz = (u(i, j, k - 2).z() - 8.0 * u(i, j, k - 1).z() +
                             8.0 * u(i, j, k + 1).z() - u(i, j, k + 2).z()) *
                           one_over_12dz;
    };
  }
} // getGradV

// Get Eddington Factor based on either constant or adaptive limiter FLD
template<std::size_t D>
void
getEddFactor(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> lambda_a,
  field<double>::accessor<wo, na> edd_factor_a,
  single<std::size_t>::accessor<ro> limiter_id_a,
  single<std::size_t>::accessor<ro> closure_id_a) noexcept {

  auto lambda = m.template mdcolex<is::cells>(lambda_a);
  auto edd_factor = m.template mdcolex<is::cells>(edd_factor_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      edd_factor(i) = AFLDEddFactor(lambda(i), *limiter_id_a, *closure_id_a);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(lambda,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      edd_factor(i, j) =
        AFLDEddFactor(lambda(i, j), *limiter_id_a, *closure_id_a);
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(lambda,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      edd_factor(i, j, k) =
        AFLDEddFactor(lambda(i, j, k), *limiter_id_a, *closure_id_a);
    };
  }

} // getEddFactor

// Get the radiation pressure tensor P
// Modified so that it takes the eddington factor field
template<std::size_t D>
void
getTensorP(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<wo,
    na> P_tensor_a,
  field<double>::accessor<ro, na> Esf_a,
  typename field<vec<D>>::template accessor<ro, na> gradEsf_a,
  field<double>::accessor<ro, na> gradE_mag_a,
  field<double>::accessor<ro, na> edd_factor_a) noexcept {

  auto P_tensor = m.template mdcolex<is::cells>(P_tensor_a);
  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto grad_esf = m.template mdcolex<is::cells>(gradEsf_a);
  auto grad_e_mag = m.template mdcolex<is::cells>(gradE_mag_a);
  auto edd_factor = m.template mdcolex<is::cells>(edd_factor_a);

  const double zero_guard = 1.0e-15;

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double f = edd_factor(i);
      P_tensor(i).xx = (0.5 * (1 - f) + 0.5 * (3 * f - 1)) * Esf(i);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(Esf,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      const double nx = grad_esf(i, j).x() / (grad_e_mag(i, j) + zero_guard);
      const double ny = grad_esf(i, j).y() / (grad_e_mag(i, j) + zero_guard);

      const double f = edd_factor(i, j);

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

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      const double nx =
        grad_esf(i, j, k).x() / (grad_e_mag(i, j, k) + zero_guard);
      const double ny =
        grad_esf(i, j, k).y() / (grad_e_mag(i, j, k) + zero_guard);
      const double nz =
        grad_esf(i, j, k).z() / (grad_e_mag(i, j, k) + zero_guard);

      const double f = edd_factor(i, j, k);

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
getGradE(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, ro> Esf_a,
  typename field<vec<D>>::template accessor<wo, ro> gradEsf_a) noexcept {
  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto grad_esf = m.template mdcolex<is::cells>(gradEsf_a);

  if constexpr(D == 1) {
    const double one_over_12dx = 1.0 / (12.0 * m.template delta<ax::x>());

    // Application of the 5-stencil central differencing:
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      grad_esf(i).x() =
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

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      // Application of the 5-stencil central differencing:
      grad_esf(i, j).x() = (Esf(i - 2, j) - 8.0 * Esf(i - 1, j) +
                             8.0 * Esf(i + 1, j) - Esf(i + 2, j)) *
                           one_over_12dx;
      grad_esf(i, j).y() = (Esf(i, j - 2) - 8.0 * Esf(i, j - 1) +
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

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      // Application of the 5-stencil central differencing:
      grad_esf(i, j, k).x() = (Esf(i - 2, j, k) - 8.0 * Esf(i - 1, j, k) +
                                8.0 * Esf(i + 1, j, k) - Esf(i + 2, j, k)) *
                              one_over_12dx;
      grad_esf(i, j, k).y() = (Esf(i, j - 2, k) - 8.0 * Esf(i, j - 1, k) +
                                8.0 * Esf(i, j + 1, k) - Esf(i, j + 2, k)) *
                              one_over_12dy;
      grad_esf(i, j, k).z() = (Esf(i, j, k - 2) - 8.0 * Esf(i, j, k - 1) +
                                8.0 * Esf(i, j, k + 1) - Esf(i, j, k + 2)) *
                              one_over_12dz;
    };
  }
} // getGradE

//
// Compute
//  1) magnitude of grad(E)
//  2) variable R (Eq 16 of Moens 2022 paper)
//  3) flux limiter function `lambda` - can either be a constant or adaptive
//
// Note : Lambda is only computed on the main grid (dm::quantities), then the
// outermost values are copied into ghost zones.
//
template<std::size_t D>
void
getLambda(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  field<double>::accessor<ro, na> Esf_a,
  typename field<vec<D>>::template accessor<ro, na> gradEsf_a,
  field<double>::accessor<rw, na> gradE_mag_a,
  field<double>::accessor<wo, na> R_a,
  field<double>::accessor<wo, na> lambda_a,
  single<double>::accessor<ro> kappa_a,
  single<std::size_t>::accessor<ro> limiter_id_a) noexcept {

  auto r = m.template mdcolex<is::cells>(r_a);
  auto Esf = m.template mdcolex<is::cells>(Esf_a);
  auto grad_esf = m.template mdcolex<is::cells>(gradEsf_a);
  auto grad_e_mag = m.template mdcolex<is::cells>(gradE_mag_a);
  auto R = m.template mdcolex<is::cells>(R_a);
  auto lambda = m.template mdcolex<is::cells>(lambda_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      auto const kappa = *kappa_a;
      const double zero_guard = 1.0e-15;

      grad_e_mag(i) = std::abs(grad_esf(i).x());
      R(i) = grad_e_mag(i) / (kappa * r(i) * Esf(i) + zero_guard);
      lambda(i) = AFLDlambda(R(i), *limiter_id_a);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(Esf,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(ji, mdpolicy_qq) {
      auto const kappa = *kappa_a;
      const double zero_guard = 1.0e-15;

      auto [j, i] = ji;
      grad_e_mag(i, j) = grad_esf(i, j).norm();
      R(i, j) = grad_e_mag(i, j) / (kappa * r(i, j) * Esf(i, j) + zero_guard);
      lambda(i, j) = AFLDlambda(R(i, j), *limiter_id_a);
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(Esf,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(kji, mdpolicy_qqq) {
      auto const kappa = *kappa_a;
      const double zero_guard = 1.0e-15;

      auto [k, j, i] = kji;
      grad_e_mag(i, j, k) = grad_esf(i, j, k).norm();
      R(i, j, k) =
        grad_e_mag(i, j, k) / (kappa * r(i, j, k) * Esf(i, j, k) + zero_guard);
      lambda(i, j, k) = AFLDlambda(R(i, j, k), *limiter_id_a);
    };
  }
} // getLambda

// Get the radiation force using the FLD approximation
template<std::size_t D>
void
getRadForce(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> lambda_a,
  typename field<vec<D>>::template accessor<ro, na> gradEsf_a,
  typename field<vec<D>>::template accessor<wo, na> fr_a) noexcept {

  auto lambda = m.template mdcolex<is::cells>(lambda_a);
  auto grad_esf = m.template mdcolex<is::cells>(gradEsf_a);
  auto fr = m.template mdcolex<is::cells>(fr_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      fr(i).x() = -lambda(i) * grad_esf(i).x();
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(lambda,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      fr(i, j).x() = -lambda(i, j) * grad_esf(i, j).x();
      fr(i, j).y() = -lambda(i, j) * grad_esf(i, j).y();
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(lambda,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());
    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      fr(i, j, k).x() = -lambda(i, j, k) * grad_esf(i, j, k).x();
      fr(i, j, k).y() = -lambda(i, j, k) * grad_esf(i, j, k).y();
      fr(i, j, k).z() = -lambda(i, j, k) * grad_esf(i, j, k).z();
    };
  }
} // getRadForce

// Radiation force and work terms. (Explicit source terms) (See Eq(21) in
// Moens2022)
template<std::size_t D>
void
explicit_source_update(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  // Primitive variables
  typename field<vec<D>>::template accessor<ro, na> velocity_a,
  // Terms required for computing radiation force and photon tiring
  typename field<vec<D>>::template accessor<ro, na> fr_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<ro,
    na> P_tensor_a,
  typename field<spec::tensor<D, spec::tensor_rank::Two>>::template accessor<ro,
    na> gradV_a,
  // time derivative
  field<double>::accessor<rw, na> dt_total_energy_density_a,
  typename field<vec<D>>::template accessor<rw, na>
    dt_momentum_energy_density_a,
  field<double>::accessor<rw, na> dt_radiation_energy_density_a) noexcept {

  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto fr = m.template mdcolex<is::cells>(fr_a);
  auto P_tensor = m.template mdcolex<is::cells>(P_tensor_a);
  auto grad_v = m.template mdcolex<is::cells>(gradV_a);

  auto dt_momentum_density =
    m.template mdcolex<is::cells>(dt_momentum_energy_density_a);
  auto dt_total_energy_density =
    m.template mdcolex<is::cells>(dt_total_energy_density_a);
  auto dt_radiation_energy_density =
    m.template mdcolex<is::cells>(dt_radiation_energy_density_a);

  // const double radiation_constant = hard::constants::cgs::radiation_constant;

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {

      // Explicitly updating conserved variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i) += fr(i);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x() * fr(i).x()
      dt_total_energy_density(i) += velocity(i).x() * fr(i).x();

      // Subtracting the photon tiring term, (P::grad_v), from the radiation
      // energy density in each cell. See Eq(34) in Moens2022.
      dt_radiation_energy_density(i) += -P_tensor(i).xx * grad_v(i).xx;

      // TODO:
      // Add the source from the temperature
      // NOTE: Isn't this already in the rad_root part?
      // dt_radiation_energy_density(i) +=
      //   constants::cgs::radiation_constant * pow(T_source(i), 4);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      // Explicitly updating conserved variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j) += fr(i, j);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x() * fr(i).x()
      dt_total_energy_density(i, j) +=
        velocity(i, j).x() * fr(i, j).x() + velocity(i, j).y() * fr(i, j).y();

      // Subtracting the photon tiring term, (P::grad_v), from the radiation
      // energy density in each cell. See Eq(34) in Moens2022.
      dt_radiation_energy_density(i, j) +=
        -(P_tensor(i, j).xx * grad_v(i, j).xx +
          P_tensor(i, j).xy * grad_v(i, j).xy +
          P_tensor(i, j).yx * grad_v(i, j).yx +
          P_tensor(i, j).yy * grad_v(i, j).yy);
    };
  }
  else {
    auto mdpolicy_qqq = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      // Explicitly updating conserved variables with S_ex
      // Adding the radiation force term to the momentum density
      dt_momentum_density(i, j, k) += fr(i, j, k);

      // Updating the total gas energy density: Adding contribution from the
      // work done by the radiative force: vdot_fr(i) = u(i).x() * fr(i).x()
      dt_total_energy_density(i, j, k) +=
        velocity(i, j, k).x() * fr(i, j, k).x() +
        velocity(i, j, k).y() * fr(i, j, k).y() +
        velocity(i, j, k).z() * fr(i, j, k).z();

      // Subtracting the photon tiring term, (P::grad_v), from the radiation
      // energy density in each cell. See Eq(34) in Moens et al. 2022.
      dt_radiation_energy_density(i, j, k) +=
        -(P_tensor(i, j, k).xx * grad_v(i, j, k).xx +
          P_tensor(i, j, k).xy * grad_v(i, j, k).xy +
          P_tensor(i, j, k).xz * grad_v(i, j, k).xz +
          P_tensor(i, j, k).yx * grad_v(i, j, k).yx +
          P_tensor(i, j, k).yy * grad_v(i, j, k).yy +
          P_tensor(i, j, k).yz * grad_v(i, j, k).yz +
          P_tensor(i, j, k).zx * grad_v(i, j, k).zx +
          P_tensor(i, j, k).zy * grad_v(i, j, k).zy +
          P_tensor(i, j, k).zz * grad_v(i, j, k).zz);
    }; // forall
  }
} // explicitSourceUpdate

// Compute the diffusion coefficients D
template<std::size_t D>
void
getDiff(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  field<double>::accessor<ro, na> lambda_a,
  field<double>::accessor<wo, na> Diff_a,
  single<double>::accessor<ro> kappa_a) noexcept {

  auto r = m.template mdcolex<is::cells>(r_a);
  auto lambda = m.template mdcolex<is::cells>(lambda_a);
  auto Diff = m.template mdcolex<is::cells>(Diff_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      auto const kappa = *kappa_a;
      const double clight = hard::constants::cgs::speed_of_light;
      Diff(i) = clight * lambda(i) / (kappa * r(i));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_aa = get_mdiota_policy(r,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_aa) {
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

    s.executor().forall(kji, mdpolicy_aaa) {
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
diffusion_init(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> Diff_a,
  typename field<double>::template accessor<wo, na> Df_xa,
  typename field<double>::template accessor<wo, na> Df_ya,
  typename field<double>::template accessor<wo, na> Df_za) noexcept {

  auto Diff = m.template mdcolex<is::cells>(Diff_a);
  auto Df_x = m.template mdcolex<is::cells>(Df_xa);
  auto Df_y = m.template mdcolex<is::cells>(Df_ya);
  auto Df_z = m.template mdcolex<is::cells>(Df_za);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::corrector>())) {
      Df_x(i) = 2 * Diff(i) * Diff(i - 1) / (Diff(i) + Diff(i - 1));
    };
  }
  else if constexpr(D == 2) {

    auto mdpolicy_cc = get_mdiota_policy(Diff,
      m.template cells<ax::y, dm::corrector>(),
      m.template cells<ax::x, dm::corrector>());

    s.executor().forall(ji, mdpolicy_cc) {
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

    s.executor().forall(kji, mdpolicy_ccc) {
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
const_init(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<wo, na> f_a,
  double w) noexcept {

  auto f = m.template mdcolex<is::cells>(f_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      f(i) = w;
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_aa = get_mdiota_policy(f,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_aa) {
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

    s.executor().forall(kji, mdpolicy_aaa) {
      auto [k, j, i] = kji;
      f(i, j, k) = w;
    };
  } // if
} // const_init

template<std::size_t D>
void
stencil_init(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> Df_xa,
  typename field<double>::template accessor<ro, ro> Df_ya,
  typename field<double>::template accessor<ro, ro> Df_za,
  typename field<stencil<D>>::template accessor<wo, na> Ew_a,
  single<double>::accessor<ro> dt_a) noexcept {
  // TODO: Stencil, Ew ghosts can be `na` (?)

  auto Df_x = m.template mdcolex<is::cells>(Df_xa);
  auto Df_y = m.template mdcolex<is::cells>(Df_ya);
  auto Df_z = m.template mdcolex<is::cells>(Df_za);
  auto Ew = m.template mdcolex<is::cells>(Ew_a);

  if constexpr(D == 1) {
    const double dx{m.template delta<ax::x>()};

    s.executor().forall(i, (m.template cells<ax::x, dm::corrector>())) {
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

    s.executor().forall(ji, mdpolicy_qq) {
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

    s.executor().forall(kji, mdpolicy_qqq) {
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

} // namespace hard::tasks::rad

#endif // HARD_MODULE_RAD_TASKS_RAD_HH
