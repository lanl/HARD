#ifndef HARD_GFM_MODULES_SHOCK_TUBE_MM_HH
#define HARD_GFM_MODULES_SHOCK_TUBE_MM_HH

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

struct oblique_sod_mm {
  static constexpr double rL = 1.0;
  static constexpr double uL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 0.125;
  static constexpr double uR = 0.0;
  static constexpr double pR = 0.1;

  // A point on the initial oblique discontinuity.
  static constexpr double x0 = 0.5;
  static constexpr double y0 = 0.5;
  static constexpr double z0 = 0.5;

  // Initialization angles for the interface normal.
  // phi   : azimuthal angle in the x-y plane.
  // theta : elevation angle out of the x-y plane (used in 3D).
  static constexpr double phi = M_PI / 10.0; // 2D and 3D
  static constexpr double theta = M_PI / 10.0; // only 3D

  static constexpr int mL = 0;
  static constexpr int mR = 1;
}; // struct oblique_sod_mm

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

struct lee_shock {
  static constexpr double rL = 0.9525;
  static constexpr double uL = 0.0;
  static constexpr double pL = 1.0;

  static constexpr double rR = 3.810;
  static constexpr double uR = 0.0;
  static constexpr double pR = 2.0;

  static constexpr double x0 = 0.5;

  static constexpr int mL = 0;
  static constexpr int mR = 1;
}; // struct lee_shock

struct helium_bubble_mm_init {
  static constexpr double rL = 1.652; // post-shocked air
  static constexpr double uL = 113.95;
  static constexpr double pL = 1.57e5;

  static constexpr double rR = 1.200; // pre-shocked air
  static constexpr double uR = 0.0;
  static constexpr double pR = 1.0e5;

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

  vec<D> uL(T::uL), uR(T::uR);

  //-----------------------------------------------------------------------------
  // TODO: singularity eos v 1.11.1 doesn't output correct internal energy for
  //       some JWL parameters, has been fixed here:
  //       https://github.com/lanl/singularity-eos/pull/645
  //       We haven't updated to that yet, so for now we use the analytical
  //       solution, reference Lee et al. 2013 (Lee Shock test) paper:
  // https://www.sciencedirect.com/science/article/pii/S0021999113002246

  spec::config_py config(opt::config.value());
  std::string eos = config["eos_list"][0]["eos"].cast<std::string>();
  std::string eos_jwl = "jwl";

  double eR, eL; // internal energy

  if(eos == eos_jwl) { // if eos is jwl
    // manual function
    spec::config_py config(opt::config.value());
    // JWL parameters from input file
    const double A = config["eos_list"][0]["JWL_A"].cast<double>();
    const double B = config["eos_list"][0]["JWL_B"].cast<double>();
    const double R1 = config["eos_list"][0]["JWL_R1"].cast<double>();
    const double R2 = config["eos_list"][0]["JWL_R2"].cast<double>();
    const double w = config["eos_list"][0]["JWL_w"].cast<double>();
    const double rho0 = config["eos_list"][0]["JWL_rho0"].cast<double>();

    double pRef_R = A * exp(-R1 * rho0 / T::rR) + B * exp(-R2 * rho0 / T::rR);
    double eRef_R = A / (rho0 * R1) * exp(-R1 * rho0 / T::rR) +
                    B / (rho0 * R2) * exp(-R2 * rho0 / T::rR);
    eR = (T::pR - pRef_R) / (w * T::rR) + eRef_R;

    double pRef_L = A * exp(-R1 * rho0 / T::rL) + B * exp(-R2 * rho0 / T::rL);
    double eRef_L = A / (rho0 * R1) * exp(-R1 * rho0 / T::rL) +
                    B / (rho0 * R2) * exp(-R2 * rho0 / T::rL);
    eL = (T::pL - pRef_L) / (w * T::rL) + eRef_L;
  }
  else {
    // singuarity function
    double sieR = 0., sieL = 0.; // initial guess, required for ieRhoPSie
    eL = eos_mm[static_cast<std::size_t>(T::mL)].ieRhoPSie(T::rL, T::pL, sieL);
    eR = eos_mm[static_cast<std::size_t>(T::mR)].ieRhoPSie(T::rR, T::pR, sieR);
  }
  //-----------------------------------------------------------------------------

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const auto x = m.template center<ax::x>(i);

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
        total_energy_density(i) = T::rL * eL + 0.5 * T::rL * (T::uL * T::uL);

        hard::tasks::util::pack_q<D>(q_scalar(i), T::rL, uL, T::pL, eL);
        q_scalar_known(i) = q_scalar(i);
      }
      else {
        material(i) = T::mR;
        frozen_material(i) = T::mR;
        psi(i) = -std::abs(phi);

        mass_density(i) = T::rR;
        momentum_density(i).x() = T::rR * T::uR;
        total_energy_density(i) = T::rR * eR + 0.5 * T::rR * (T::uR * T::uR);

        hard::tasks::util::pack_q<D>(q_scalar(i), T::rR, uR, T::pR, eR);
        q_scalar_known(i) = q_scalar(i);

      } // if

      // Initial exact constant extrapolation fields.
      // plus: left/minus material extended into plus/right side.
      hard::tasks::util::pack_q<D>(q_ext_plus(i), T::rL, uL, T::pL, eL);

      // minus: right/plus material extended into minus/left side.
      hard::tasks::util::pack_q<D>(q_ext_minus(i), T::rR, uR, T::pR, eR);

    } // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto ji : mdpolicy_yx) {
      auto [j, i] = ji;
      const auto x = m.template center<ax::x>(i);

      const double phi = x - T::x0;

      if(x < T::x0) {
        material(i, j) = T::mL;
        frozen_material(i, j) = T::mL;
        psi(i, j) = -std::abs(phi);

        mass_density(i, j) = T::rL;
        momentum_density(i, j) = vec<D>(T::rL * T::uL);
        total_energy_density(i, j) = T::rL * eL + 0.5 * T::rL * (T::uL * T::uL);

        hard::tasks::util::pack_q<D>(q_scalar(i, j), T::rL, uL, T::pL, eL);
        q_scalar_known(i, j) = q_scalar(i, j);
      }
      else {
        material(i, j) = T::mR;
        frozen_material(i, j) = T::mR;
        psi(i, j) = -std::abs(phi);

        mass_density(i, j) = T::rR;
        momentum_density(i, j) = vec<D>(T::rR * T::uR);
        total_energy_density(i, j) = T::rR * eR + 0.5 * T::rR * (T::uR * T::uR);

        hard::tasks::util::pack_q<D>(q_scalar(i, j), T::rR, uR, T::pR, eR);
        q_scalar_known(i, j) = q_scalar(i, j);
      } // if

      // Initial exact constant extrapolation fields.
      // plus: left/minus material extended into plus/right side.
      hard::tasks::util::pack_q<D>(q_ext_plus(i, j), T::rL, uL, T::pL, eL);

      // minus: right/plus material extended into minus/left side.
      hard::tasks::util::pack_q<D>(q_ext_minus(i, j), T::rR, uR, T::pR, eR);
    } // for
  }
  else /* D == 3 */ {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto kji : mdpolicy_zyx) {
      auto [k, j, i] = kji;
      const auto x = m.template center<ax::x>(i);

      const double phi = x - T::x0;

      if(x < T::x0) {
        material(i, j, k) = T::mL;
        frozen_material(i, j, k) = T::mL;
        psi(i, j, k) = -std::abs(phi);

        mass_density(i, j, k) = T::rL;
        momentum_density(i, j, k) = vec<D>(T::rL * T::uL);
        total_energy_density(i, j, k) =
          T::rL * eL + 0.5 * T::rL * (T::uL * T::uL);
        hard::tasks::util::pack_q<D>(q_scalar(i, j, k), T::rL, uL, T::pL, eL);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      }
      else {
        material(i, j, k) = T::mR;
        frozen_material(i, j, k) = T::mR;
        psi(i, j, k) = -std::abs(phi);

        mass_density(i, j, k) = T::rR;
        momentum_density(i, j, k) = vec<D>(T::rR * T::uR);
        total_energy_density(i, j, k) =
          T::rR * eR + 0.5 * T::rR * (T::uR * T::uR);
        hard::tasks::util::pack_q<D>(q_scalar(i, j, k), T::rR, uR, T::pR, eR);
        q_scalar_known(i, j, k) = q_scalar(i, j, k);
      } // if

      // Initial exact constant extrapolation fields.
      // plus: left/minus material extended into plus/right side.
      hard::tasks::util::pack_q<D>(q_ext_plus(i, j, k), T::rL, uL, T::pL, eL);

      // minus: right/plus material extended into minus/left side.
      hard::tasks::util::pack_q<D>(q_ext_minus(i, j, k), T::rR, uR, T::pR, eR);
    } // for
  } // if
} // shock

} // namespace hard::tasks::initial_data

#endif // HARD_GFM_MODULES_SHOCK_TUBE_MM_HH
