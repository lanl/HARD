#ifndef HARD_GFM_MODULES_BUBBLE_HH
#define HARD_GFM_MODULES_BUBBLE_HH

#include "../modules/common/tasks/utils.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

/*----------------------------------------------------------------------------*
    Bubble Shock Interaction.
 *----------------------------------------------------------------------------*/

namespace bubbles {

struct he_bubble {
  static constexpr double rL = 1.337;
  static constexpr double uL = 1.22; // shocked air
  static constexpr double pL = 1.19;

  static constexpr double rM = 1.0; // unshocked air
  static constexpr double uM = 0.0;
  static constexpr double pM = 1.0;

  static constexpr double rR = 0.13798; // bubble
  static constexpr double uR = 0.0;
  static constexpr double pR = 1.0;

  static constexpr double x0 = 0.5; // position of the shock

  static constexpr int mL = 0; // air
  static constexpr int mM = 0;
  static constexpr int mR = 1; // bubble
}; // struct he bubble

struct he_bubble_extra {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 10.0;

  static constexpr double rM = 1.0;
  static constexpr double uM = 0.0;
  static constexpr double pM = 1.0;

  static constexpr double rR = 0.13798; // bubble
  static constexpr double uR = 0.0;
  static constexpr double pR = 1.0;

  static constexpr double x0 = 0.5; // position of the shock

  static constexpr int mL = 0; // air
  static constexpr int mM = 0;
  static constexpr int mR = 1; // bubble
}; // struct he bubble

} // namespace bubbles

template<typename T, std::size_t D>
auto
he_bubble_mm(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  typename field<int>::accessor<wo, na> material_a,
  typename field<int>::accessor<wo, na> frozen_material_a,
  field<double>::accessor<wo, na> psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_scalar_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_scalar_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_ext_plus_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_ext_minus_a,
  const std::vector<eos::eos_wrapper> & eos_mm) {

  using common::tasks::utils::get_mdiota_policy;

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto material = m.template mdcolex<is::cells>(material_a);
  auto frozen_material = m.template mdcolex<is::cells>(frozen_material_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto q_scalar = m.template mdcolex<is::cells>(q_scalar_a);
  auto q_scalar_known = m.template mdcolex<is::cells>(q_scalar_known_a);
  auto q_ext_plus = m.template mdcolex<is::cells>(q_ext_plus_a);
  auto q_ext_minus = m.template mdcolex<is::cells>(q_ext_minus_a);

  // initial conditions
  const double eL = util::find_sie(eos_mm[T::mL], T::rL, T::pL);
  const double eM = util::find_sie(eos_mm[T::mM], T::rM, T::pM);
  const double eR = util::find_sie(eos_mm[T::mR], T::rR, T::pR);

  if constexpr(D == 1) {
    flog_fatal("Shocked bubble problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    vec<D> uL, uM, uR;
    uL.x() = T::uL;
    uL.y() = 0.0;
    uM.x() = T::uM;
    uM.y() = 0.0;
    uR.x() = T::uR;
    uR.y() = 0.0;

    constexpr double xb = 1.5;
    constexpr double yb = 1.0;
    constexpr double rb = 0.5;

    for(auto ji : mdpolicy_yx) {
      auto [j, i] = ji;
      const auto x = m.template center<ax::x>(i);
      const auto y = m.template center<ax::y>(j);

      const double dx = x - xb;
      const double dy = y - yb;
      const double s = std::sqrt(dx * dx + dy * dy);

      // Material-signed distance used by the HARD level-set update:
      // psi must be negative on both sides, with magnitude equal to the
      // distance to the closest interface. The sign relative to a neighbor of a
      // different material is handled later by rhs_upwind via sign flipping.
      psi(i, j) = -std::abs(s - rb);

      if(s <= rb) {
        material(i, j) = T::mR;
        frozen_material(i, j) = T::mR;

        mass_density(i, j) = T::rR;
        momentum_density(i, j).x() = T::rR * T::uR;
        momentum_density(i, j).y() = 0.0;
        total_energy_density(i, j) = T::rR * eR + 0.5 * T::rR * (T::uR * T::uR);

        common::tasks::utils::pack_q<D>(q_scalar(i, j), T::rR, uR, T::pR, eR);
        q_scalar_known(i, j) = q_scalar(i, j);
      }
      else if(x < T::x0) {
        // Shocked / left state
        material(i, j) = T::mL;
        frozen_material(i, j) = T::mL;

        mass_density(i, j) = T::rL;
        momentum_density(i, j).x() = T::rL * T::uL;
        momentum_density(i, j).y() = 0.0;
        total_energy_density(i, j) = T::rL * eL + 0.5 * T::rL * (T::uL * T::uL);

        common::tasks::utils::pack_q<D>(q_scalar(i, j), T::rL, uL, T::pL, eL);
        q_scalar_known(i, j) = q_scalar(i, j);
      }
      else {
        // Unshocked ambient / middle state
        material(i, j) = T::mM;
        frozen_material(i, j) = T::mM;

        mass_density(i, j) = T::rM;
        momentum_density(i, j).x() = T::rM * T::uM;
        momentum_density(i, j).y() = 0.0;
        total_energy_density(i, j) = T::rM * eM + 0.5 * T::rM * (T::uM * T::uM);

        common::tasks::utils::pack_q<D>(q_scalar(i, j), T::rM, uM, T::pM, eM);
        q_scalar_known(i, j) = q_scalar(i, j);
      }

      // Extrapolated states across the bubble interface
      common::tasks::utils::pack_q<D>(q_ext_plus(i, j), T::rM, uM, T::pM, eM);
      common::tasks::utils::pack_q<D>(q_ext_minus(i, j), T::rR, uR, T::pR, eR);
    }
  }
  else if constexpr(D == 3) {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    vec<D> uL, uM, uR;
    uL.x() = T::uL;
    uL.y() = 0.0;
    uL.z() = 0.0;
    uM.x() = T::uM;
    uM.y() = 0.0;
    uM.z() = 0.0;
    uR.x() = T::uR;
    uR.y() = 0.0;
    uR.z() = 0.0;

    constexpr double xb = 1.5;
    constexpr double yb = 1.0;
    constexpr double zb = yb;
    constexpr double rb = 0.5;

    for(auto kji : mdpolicy_zyx) {
      auto [k, j, i] = kji;
      const auto x = m.template center<ax::x>(i);
      const auto y = m.template center<ax::y>(j);
      const auto z = m.template center<ax::z>(k);

      const double dx = x - xb;
      const double dy = y - yb;
      const double dz = z - zb;
      const double s = std::sqrt(dx * dx + dy * dy + dz * dz);

      // Material-signed distance used by the HARD level-set update:
      // psi must be negative on both sides, with magnitude equal to the
      // distance to the closest interface. The sign relative to a neighbor of
      // a different material is handled later by rhs_upwind via sign flipping.
      psi(i, j, k) = -std::abs(s - rb);

      if(s <= rb) {
        material(i, j, k) = T::mR;
        frozen_material(i, j, k) = T::mR;

        mass_density(i, j, k) = T::rR;
        momentum_density(i, j, k).x() = T::rR * T::uR;
        momentum_density(i, j, k).y() = 0.0;
        momentum_density(i, j, k).z() = 0.0;
        total_energy_density(i, j, k) =
          T::rR * eR + 0.5 * T::rR * (T::uR * T::uR);

        common::tasks::utils::pack_q<D>(
          q_scalar(i, j, k), T::rR, uR, T::pR, eR);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      }
      else if(x < T::x0) {
        // Shocked / left state
        material(i, j, k) = T::mL;
        frozen_material(i, j, k) = T::mL;

        mass_density(i, j, k) = T::rL;
        momentum_density(i, j, k).x() = T::rL * T::uL;
        momentum_density(i, j, k).y() = 0.0;
        momentum_density(i, j, k).z() = 0.0;
        total_energy_density(i, j, k) =
          T::rL * eL + 0.5 * T::rL * (T::uL * T::uL);

        common::tasks::utils::pack_q<D>(
          q_scalar(i, j, k), T::rL, uL, T::pL, eL);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      }
      else {
        // Unshocked ambient / middle state
        material(i, j, k) = T::mM;
        frozen_material(i, j, k) = T::mM;

        mass_density(i, j, k) = T::rM;
        momentum_density(i, j, k).x() = T::rM * T::uM;
        momentum_density(i, j, k).y() = 0.0;
        momentum_density(i, j, k).z() = 0.0;
        total_energy_density(i, j, k) =
          T::rM * eM + 0.5 * T::rM * (T::uM * T::uM);

        common::tasks::utils::pack_q<D>(
          q_scalar(i, j, k), T::rM, uM, T::pM, eM);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      }

      // Extrapolated states across the bubble interface
      common::tasks::utils::pack_q<D>(
        q_ext_plus(i, j, k), T::rM, uM, T::pM, eM);
      common::tasks::utils::pack_q<D>(
        q_ext_minus(i, j, k), T::rR, uR, T::pR, eR);
    }
  }

} // he_bubble_mm

} // namespace hard::tasks::initial_data

#endif // HARD_GFM_MODULES_BUBBLE_HH