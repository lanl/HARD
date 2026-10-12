#ifndef HARD_GFM_MODULES_OBLIQUE_SOD_MM_HH
#define HARD_GFM_MODULES_OBLIQUE_SOD_MM_HH

#include "../modules/common/tasks/utils.hh"
#include "../utils.hh"
#include "shock_tube_mm.hh"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace hard::tasks::initial_data {

struct oblique_shock_mm_task : annotation_task_region {
  inline static const std::string name{"oblique_shock_mm"};
};

template<typename T, std::size_t D>
auto
oblique_shock_mm(flecsi::exec::cpu s,
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
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {

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

  spec::config_py config(opt::config.value());
  std::string eos = config["eos_list"][0]["eos"].cast<std::string>();
  std::string eos_jwl = "jwl";

  double eR, eL;

  if(eos == eos_jwl) {
    spec::config_py config(opt::config.value());
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
    double sieR = 0., sieL = 0.;
    eL = eos_mm[static_cast<std::size_t>(T::mL)].ieRhoPSie(T::rL, T::pL, sieL);
    eR = eos_mm[static_cast<std::size_t>(T::mR)].ieRhoPSie(T::rR, T::pR, sieR);
  }

  const double nx = std::cos(T::theta) * std::cos(T::phi);
  const double ny = std::cos(T::theta) * std::sin(T::phi);
  const double nz = std::sin(T::theta);

  if constexpr(D == 1) {
    flog_fatal("oblique_sod_mm only runs in 2D and 3D; for the "
               "1D case use normal sod_mm");
  }
  else if constexpr(D == 2) {
    auto mdpolicy_yx = get_mdiota_policy(mass_density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto ji : mdpolicy_yx) {
      auto [j, i] = ji;
      const auto x = m.template center<ax::x>(i);
      const auto y = m.template center<ax::y>(j);
      const double phi = nx * (x - T::x0) + ny * (y - T::y0);

      if(phi < 0.0) {
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
      }

      hard::tasks::util::pack_q<D>(q_ext_plus(i, j), T::rL, uL, T::pL, eL);
      hard::tasks::util::pack_q<D>(q_ext_minus(i, j), T::rR, uR, T::pR, eR);
    }
  }
  else {
    auto mdpolicy_zyx = get_mdiota_policy(mass_density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    for(auto kji : mdpolicy_zyx) {
      auto [k, j, i] = kji;
      const auto x = m.template center<ax::x>(i);
      const auto y = m.template center<ax::y>(j);
      const auto z = m.template center<ax::z>(k);
      const double phi = nx * (x - T::x0) + ny * (y - T::y0) + nz * (z - T::z0);

      if(phi < 0.0) {
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
      }

      hard::tasks::util::pack_q<D>(q_ext_plus(i, j, k), T::rL, uL, T::pL, eL);
      hard::tasks::util::pack_q<D>(q_ext_minus(i, j, k), T::rR, uR, T::pR, eR);
    }
  }
}

} // namespace hard::tasks::initial_data

#endif // HARD_GFM_MODULES_OBLIQUE_SOD_MM_HH
