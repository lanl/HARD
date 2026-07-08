#ifndef HARD_GFM_MODULES_SHOCK_TUBE_MM_HH
#define HARD_GFM_MODULES_SHOCK_TUBE_MM_HH

#include "../modules/common/tasks/utils.hh"
#include "../utils.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

namespace shock_tubes {

struct sod_mm {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 0.0;
  static constexpr double pR = 0.1;

  static constexpr double x0 = 0.5; // 0.5

  static constexpr int mL = 0;
  static constexpr int mR = 1;
}; // struct sod_mm

struct shyue_shock {
  static constexpr double rL = 1.7;
  static constexpr double uL = 0.0;
  static constexpr double pL = 10.0;

  static constexpr double rR = 1.0;
  static constexpr double uR = 0.0;
  static constexpr double pR = 0.5;

  static constexpr double x0 = 0.5;

  static constexpr int mL = 0;
  static constexpr int mR = 1;
}; // struct shyue_shock

struct helium_bubble_mm_init {
  static constexpr double rL = 1.7;
  static constexpr double uL = 0.0;
  static constexpr double pL = 10.0;

  static constexpr double rR = 1.0;
  static constexpr double uR = 0.0;
  static constexpr double pR = 0.5;

  static constexpr double x0 = 100.0;

  // Both sides of the initial discontinuity are ambient gas.
  static constexpr int mL = 0;
  static constexpr int mR = 0;
}; // struct helium_bubble_mm

struct CONTACT_stat {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 0.0;
  static constexpr double pR = 1.0;

  static constexpr double x0 = 0.5;

  static constexpr int mL = 0;
  static constexpr int mR = 1;
}; // struct stationary contact test

struct CONTACT_move {
  static constexpr double rL = 1.0;
  static constexpr double uL = 1.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 1.0;
  static constexpr double pR = 1.0;

  static constexpr double x0 = 0.5;

  static constexpr int mL = 0;
  static constexpr int mR = 1;
}; // struct moving contact test

} // namespace shock_tubes

/*----------------------------------------------------------------------------*
  Shock Tube.
 *----------------------------------------------------------------------------*/

struct shock_mm_task : annotation_task_region {
  inline static const std::string name{"shock_mm"};
};

template<typename T, std::size_t D>
auto
shock_mm(flecsi::exec::cpu s,
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

  auto pack_q = [](auto & qcell,
                  double rho,
                  const vec<D> & velocity,
                  double pressure,
                  double sie) {
    qcell.rho() = rho;
    qcell.u() = velocity;
    qcell.p() = pressure;
    qcell.e() = sie;
    qcell.n_comps() = 0;
  };

  vec<D> uL(T::uL), uR(T::uR);
  double sieR, sieL, sie;
  const double eL =
    eos_mm[static_cast<std::size_t>(T::mL)].ieRhoPSie(T::rL, T::pL, sieL);
  const double eR =
    eos_mm[static_cast<std::size_t>(T::mR)].ieRhoPSie(T::rR, T::pR, sieR);

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const auto x = m.template right<ax::x>(i);

      // Conventional signed level set for the shock-tube interface.
      // phi < 0 : left material
      // phi > 0 : right material
      const double phi = x - T::x0;

      if(phi < 0.0) {
        material(i) = T::mL;
        frozen_material(i) = T::mL;
        psi(i) = -std::abs(phi);

        mass_density(i) = T::rL;
        momentum_density(i).x() = T::rL * T::uL;
        const double e =
          eos_mm[static_cast<std::size_t>(T::mL)].ieRhoPSie(T::rL, T::pL, sieL);
        total_energy_density(i) = T::rL * e + 0.5 * T::rL * (T::uL * T::uL);

        pack_q(q_scalar(i), T::rL, uL, T::pL, eL);
        q_scalar_known(i) = q_scalar(i);
      }
      else {
        material(i) = T::mR;
        frozen_material(i) = T::mR;
        psi(i) = -std::abs(phi);

        mass_density(i) = T::rR;
        momentum_density(i).x() = T::rR * T::uR;
        const double e =
          eos_mm[static_cast<std::size_t>(T::mR)].ieRhoPSie(T::rR, T::pR, sieR);
        total_energy_density(i) = T::rR * e + 0.5 * T::rR * (T::uR * T::uR);

        pack_q(q_scalar(i), T::rR, uR, T::pR, eR);
        q_scalar_known(i) = q_scalar(i);

      } // if

      // Initial exact constant extrapolation fields.
      // plus: left/minus material extended into plus/right side.
      pack_q(q_ext_plus(i), T::rL, uL, T::pL, eL);

      // minus: right/plus material extended into minus/left side.
      pack_q(q_ext_minus(i), T::rR, uR, T::pR, eR);

    } // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto ji : mdpolicy_yx) {
      auto [j, i] = ji;
      const auto x = m.template right<ax::x>(i);

      const double phi = x - T::x0;

      if(x < T::x0) {
        material(i, j) = T::mL;
        frozen_material(i, j) = T::mL;
        psi(i, j) = -std::abs(phi);

        mass_density(i, j) = T::rL;
        momentum_density(i, j) = vec<D>(T::rL * T::uL);
        const double e =
          eos_mm[static_cast<std::size_t>(T::mL)].ieRhoPSie(T::rL, T::pL, sieL);
        total_energy_density(i, j) = T::rL * e + 0.5 * T::rL * (T::uL * T::uL);

        pack_q(q_scalar(i, j), T::rL, uL, T::pL, eL);
        q_scalar_known(i, j) = q_scalar(i, j);
      }
      else {
        material(i, j) = T::mR;
        frozen_material(i, j) = T::mR;
        psi(i, j) = -std::abs(phi);

        mass_density(i, j) = T::rR;
        momentum_density(i, j) = vec<D>(T::rR * T::uR);
        const double e =
          eos_mm[static_cast<std::size_t>(T::mR)].ieRhoPSie(T::rR, T::pR, sieR);
        total_energy_density(i, j) = T::rR * e + 0.5 * T::rR * (T::uR * T::uR);

        pack_q(q_scalar(i, j), T::rR, uR, T::pR, eR);
        q_scalar_known(i, j) = q_scalar(i, j);
      } // if

      // Initial exact constant extrapolation fields.
      // plus: left/minus material extended into plus/right side.
      pack_q(q_ext_plus(i, j), T::rL, uL, T::pL, eL);

      // minus: right/plus material extended into minus/left side.
      pack_q(q_ext_minus(i, j), T::rR, uR, T::pR, eR);
    } // for
  }
  else /* D == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto kji : mdpolicy_zyx) {
      auto [k, j, i] = kji;
      const auto x = m.template right<ax::x>(i);

      const double phi = x - T::x0;

      if(x < T::x0) {
        material(i, j, k) = T::mL;
        frozen_material(i, j, k) = T::mL;
        psi(i, j, k) = -std::abs(phi);

        mass_density(i, j, k) = T::rL;
        momentum_density(i, j, k) = vec<D>(T::rL * T::uL);
        const double e =
          eos_mm[static_cast<std::size_t>(T::mL)].ieRhoPSie(T::rL, T::pL, sieL);
        total_energy_density(i, j, k) =
          T::rL * e + 0.5 * T::rL * (T::uL * T::uL);
        pack_q(q_scalar(i, j, k), T::rL, uL, T::pL, eL);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      }
      else {
        material(i, j, k) = T::mR;
        frozen_material(i, j, k) = T::mR;
        psi(i, j, k) = -std::abs(phi);

        mass_density(i, j, k) = T::rR;
        momentum_density(i, j, k) = vec<D>(T::rR * T::uR);
        const double e =
          eos_mm[static_cast<std::size_t>(T::mR)].ieRhoPSie(T::rR, T::pR, sieR);
        total_energy_density(i, j, k) =
          T::rR * e + 0.5 * T::rR * (T::uR * T::uR);
        pack_q(q_scalar(i, j, k), T::rR, uR, T::pR, eR);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      } // if

      // Initial exact constant extrapolation fields.
      // plus: left/minus material extended into plus/right side.
      pack_q(q_ext_plus(i, j, k), T::rL, uL, T::pL, eL);

      // minus: right/plus material extended into minus/left side.
      pack_q(q_ext_minus(i, j, k), T::rR, uR, T::pR, eR);
    } // for
  } // if
} // shock

} // namespace hard::tasks::initial_data

#endif // HARD_GFM_MODULES_SHOCK_TUBE_MM_HH
