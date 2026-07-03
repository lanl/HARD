#ifndef HARD_MODULES_GFM_TASKS_LS_CORE_HH
#define HARD_MODULES_GFM_TASKS_LS_CORE_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/mesh.hh"
#include "../state.hh"

using namespace hard::tasks::util;

namespace hard::tasks::level_set {

struct freeze_material_task : annotation_task_region {
  inline static const std::string name{"freeze_material"};
};

// Freeze material labels at the start of each step
template<std::size_t D>
void
freeze_material(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> material_a,
  typename field<int>::template accessor<wo, ro> frozen_material_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  s.executor().forall(
    i, flecsi::util::iota_view({}, frozen_material_a.span().size())) {
    frozen_material_a(i) = material_a(i);
  };
}

struct rhs_upwind_task : annotation_task_region {
  inline static const std::string name{"rhs_upwind"};
};

// upwind HJ RHS (acts on material-signed psi_m)
template<std::size_t D>
void
rhs_upwind(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> velocity_a,
  typename field<int>::template accessor<ro, ro> frozen_material_a,
  typename field<double>::template accessor<wo, ro> rhs_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto vel = m.template mdcolex<is::cells>(velocity_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const int ip1 = i + 1;
      const int im1 = i - 1;

      auto [um, up] = minmax(0.0, vel(i));

      // material-signed neighbor values relative to center frozen material
      const double pm_ip1 = (fmat(ip1) == fmat(i)) ? psi(ip1) : -psi(ip1);
      const double pm_im1 = (fmat(im1) == fmat(i)) ? psi(im1) : -psi(im1);
      const double pm_i = psi(i);

      // adv_x in flux-form upwind; rhs = -adv to match working HJ_upwind_L_2D
      const double adv_x =
        (pm_i * (up.x() - um.x()) + pm_ip1 * um.x() - pm_im1 * up.x()) /
        m.template delta<ax::x>();

      rhs(i) = -adv_x;
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      auto [um, up] = minmax(0.0, vel(i, j));

      // Compute material-signed neighbors
      std::array<std::tuple<int, int>, 2 * D> neighbor_coords = {
        {{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}}};

      std::array<double, 2 * D> pm;
      for(int n = 0; n < 2 * D; ++n) {
        auto [ni, nj] = neighbor_coords[n];
        pm[n] = (fmat(ni, nj) == fmat(i, j)) ? psi(ni, nj) : -psi(ni, nj);
      }
      const double pm_ij = psi(i, j);

      // pm[0]=im1, pm[1]=ip1, pm[2]=jm1, pm[3]=jp1
      const double adv_x =
        (pm_ij * (up.x() - um.x()) + pm[1] * um.x() - pm[0] * up.x()) /
        m.template delta<ax::x>();

      const double adv_y =
        (pm_ij * (up.y() - um.y()) + pm[3] * um.y() - pm[2] * up.y()) /
        m.template delta<ax::y>();

      rhs(i, j) = -(adv_x + adv_y);
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(vel,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      auto [um, up] = minmax(0.0, vel(i, j, k));

      // Compute material-signed neighbors
      std::array<std::tuple<int, int, int>, 2 * D> neighbor_coords = {
        {{i - 1, j, k},
          {i + 1, j, k},
          {i, j - 1, k},
          {i, j + 1, k},
          {i, j, k - 1},
          {i, j, k + 1}}};

      std::array<double, 2 * D> pm;
      for(int n = 0; n < 2 * D; ++n) {
        auto [ni, nj, nk] = neighbor_coords[n];
        pm[n] = (fmat(ni, nj, nk) == fmat(i, j, k)) ? psi(ni, nj, nk)
                                                    : -psi(ni, nj, nk);
      }
      const double pm_ijk = psi(i, j, k);

      // pm[0]=im1, pm[1]=ip1, pm[2]=jm1, pm[3]=jp1, pm[4]=km1, pm[5]=kp1
      const double adv_x =
        (pm_ijk * (up.x() - um.x()) + pm[1] * um.x() - pm[0] * up.x()) /
        m.template delta<ax::x>();
      const double adv_y =
        (pm_ijk * (up.y() - um.y()) + pm[3] * um.y() - pm[2] * up.y()) /
        m.template delta<ax::y>();
      const double adv_z =
        (pm_ijk * (up.z() - um.z()) + pm[5] * um.z() - pm[4] * up.z()) /
        m.template delta<ax::z>();

      rhs(i, j, k) = -(adv_x + adv_y + adv_z);
    };
  }
}

struct advance_psi_s1_task : annotation_task_region {
  inline static const std::string name{"advance_psi_s1"};
};

// SSPRK3 stage update
template<std::size_t D, std::size_t STEP>
void
advance_psi_s1(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_n0_a,
  typename field<double>::template accessor<ro, ro> rhs_a,
  typename field<double>::template accessor<wo, ro> psi_out_a,
  typename single<double>::template accessor<ro> dt_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto psi_n0 = m.template mdcolex<is::cells>(psi_n0_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);
  auto psi_out = m.template mdcolex<is::cells>(psi_out_a);

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      psi_out(i) = psi_n0(i) + *dt_a * rhs(i);
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      psi_out(i, j) = psi_n0(i, j) + *dt_a * rhs(i, j);
    };
  }
  else {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      psi_out(i, j, k) = psi_n0(i, j, k) + *dt_a * rhs(i, j, k);
    };
  }
}

struct advance_psi_task : annotation_task_region {
  inline static const std::string name{"advance_psi"};
};

template<std::size_t D, std::size_t STEP>
void
advance_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> psi_n0_a,
  typename field<double>::template accessor<ro, ro> rhs_a,
  typename field<double>::template accessor<ro, ro> psi_prev_a,
  typename field<double>::template accessor<wo, ro> psi_out_a,
  typename single<double>::template accessor<ro> dt_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  static_assert(STEP >= 2 && STEP <= 3, "STEP must be 2 or 3 here");

  auto psi_n0 = m.template mdcolex<is::cells>(psi_n0_a);
  auto rhs = m.template mdcolex<is::cells>(rhs_a);
  auto psi_prev = m.template mdcolex<is::cells>(psi_prev_a);
  auto psi_out = m.template mdcolex<is::cells>(psi_out_a);

  // const double dt = *dt_a;

  auto update_one = [=](double psi0, double psip, double k) {
    if constexpr(STEP == 2) {
      return 0.75 * psi0 + 0.25 * (psip + *dt_a * k);
    }
    else {
      return (1.0 / 3.0) * psi0 + (2.0 / 3.0) * (psip + *dt_a * k);
    }
  };

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      psi_out(i) = update_one(psi_n0(i), psi_prev(i), rhs(i));
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      psi_out(i, j) = update_one(psi_n0(i, j), psi_prev(i, j), rhs(i, j));
    };
  }
  else {
    auto pol = get_mdiota_policy(psi_out,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      psi_out(i, j, k) =
        update_one(psi_n0(i, j, k), psi_prev(i, j, k), rhs(i, j, k));
    };
  }
}

struct determine_psi_task : annotation_task_region {
  inline static const std::string name{"determine_psi"};
};

// Determine state at time level n+1 (2-neigh 1D, 4-neigh 2D, 6-neigh 3D)
// algorithm 1 in Aslam et al. 2024
template<std::size_t D>
void
determine_psi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<ro, ro> frozen_psi_a, // psi3
  typename field<double>::template accessor<wo, ro> psi_a, // psi^{n+1}
  typename field<int>::template accessor<ro, ro> frozen_material_a, // L^n
  typename field<int>::template accessor<wo, ro> material_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto fpsi = m.template mdcolex<is::cells>(frozen_psi_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto fmat = m.template mdcolex<is::cells>(frozen_material_a);
  auto mat = m.template mdcolex<is::cells>(material_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      if(fpsi(i) > 0.0) {
        double psimin = fpsi(i);
        int lmin = fmat(i);

        // Array of neighbor indices
        std::array<unsigned int, 2 * D> neighbors = {i - 1, i + 1};

        for(auto ni : neighbors) {
          if(fmat(ni) != fmat(i) && fpsi(ni) < psimin) {
            psimin = fpsi(ni);
            lmin = fmat(ni);
          }
        }

        psi(i) = (psimin >= 0.0 ? 1.0 : -1.0) * fpsi(i);
        mat(i) = lmin;
      }
      else {
        psi(i) = fpsi(i);
        mat(i) = fmat(i);
      }
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      if(fpsi(i, j) > 0.0) {
        double psimin = fpsi(i, j);
        int lmin = fmat(i, j);

        // Array of neighbor coordinate pairs
        std::array<std::tuple<unsigned int, unsigned int>, 2 * D> neighbors = {
          {{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}}};

        for(auto [ni, nj] : neighbors) {
          if(fmat(ni, nj) != fmat(i, j) && fpsi(ni, nj) < psimin) {
            psimin = fpsi(ni, nj);
            lmin = fmat(ni, nj);
          }
        }

        psi(i, j) = (psimin >= 0.0 ? 1.0 : -1.0) * fpsi(i, j);
        mat(i, j) = lmin;
      }
      else {
        psi(i, j) = fpsi(i, j);
        mat(i, j) = fmat(i, j);
      }
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      if(fpsi(i, j, k) > 0.0) {
        double psimin = fpsi(i, j, k);
        int lmin = fmat(i, j, k);

        // Array of neighbor coordinate triples
        std::array<std::tuple<unsigned int, unsigned int, unsigned int>, 2 * D>
          neighbors = {{{i - 1, j, k},
            {i + 1, j, k},
            {i, j - 1, k},
            {i, j + 1, k},
            {i, j, k - 1},
            {i, j, k + 1}}};

        for(auto [ni, nj, nk] : neighbors) {
          if(fmat(ni, nj, nk) != fmat(i, j, k) && fpsi(ni, nj, nk) < psimin) {
            psimin = fpsi(ni, nj, nk);
            lmin = fmat(ni, nj, nk);
          }
        }

        psi(i, j, k) = (psimin >= 0.0 ? 1.0 : -1.0) * fpsi(i, j, k);
        mat(i, j, k) = lmin;
      }
      else {
        psi(i, j, k) = fpsi(i, j, k);
        mat(i, j, k) = fmat(i, j, k);
      }
    };
  }
}

template<std::size_t D>
FLECSI_INLINE_TARGET double
cell_sie_from_cons(double rho, const vec<D> & mom, double rhoE) noexcept {
  const vec<D> u = (1.0 / rho) * mom;
  return (rhoE - 0.5 * rho * u.norm_squared()) / rho;
}

FLECSI_INLINE_TARGET bool
rel_close_local(double a,
  double b,
  double rtol = 1.0e-10,
  double atol = 1.0e-12) noexcept {
  return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

template<std::size_t D, typename eosType>
FLECSI_INLINE_TARGET bool
eos_pair_equivalent_at_cell(const eosType & eos_disp,
  int old_mat,
  int new_mat,
  double rho,
  const vec<D> & mom,
  double rhoE) noexcept {

  if(old_mat == new_mat) {
    return true;
  }

  if(old_mat < 0 || new_mat < 0 ||
     static_cast<std::size_t>(old_mat) >= eos_disp.size() ||
     static_cast<std::size_t>(new_mat) >= eos_disp.size()) {
    return false;
  }

  if(!std::isfinite(rho) || rho <= 0.0) {
    return false;
  }

  const double e = cell_sie_from_cons<D>(rho, mom, rhoE);

  if(!std::isfinite(e) || e <= 0.0) {
    return false;
  }

  const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
  const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

  const double p_old = eos_old.pRhoSie(rho, e);
  const double p_new = eos_new.pRhoSie(rho, e);
  const double c_old = eos_old.cRhoSie(rho, e);
  const double c_new = eos_new.cRhoSie(rho, e);

  return rel_close_local(p_old, p_new) && rel_close_local(c_old, c_new);
}

struct promote_changed_material_cells_task : annotation_task_region {
  inline static const std::string name{"promote_changed_material_cells"};
};

template<std::size_t D>
void
promote_changed_material_cells(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> old_material_a,
  typename field<int>::template accessor<ro, ro> new_material_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro>
    q_ext_plus_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, ro>
    q_ext_minus_a,
  field<double>::accessor<rw, na> mass_density_a,
  typename field<vec<D>>::template accessor<rw, na> momentum_density_a,
  field<double>::accessor<rw, na> total_energy_density_a,
  typename field<vec<D>>::template accessor<rw, na> velocity_a,
  field<double>::accessor<rw, na> pressure_a,
  field<double>::accessor<rw, na> sie_a,
  field<double>::accessor<rw, na> sound_speed_a,
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {

  eos::eos_dispatch<flecsi::exec::accelerator> eos_disp(eos_mm);

  using common::tasks::utils::get_mdiota_policy;

  auto old_material = m.template mdcolex<is::cells>(old_material_a);
  auto new_material = m.template mdcolex<is::cells>(new_material_a);

  auto q_ext_plus = m.template mdcolex<is::cells>(q_ext_plus_a);
  auto q_ext_minus = m.template mdcolex<is::cells>(q_ext_minus_a);

  auto rho = m.template mdcolex<is::cells>(mass_density_a);
  auto mom = m.template mdcolex<is::cells>(momentum_density_a);
  auto rhoE = m.template mdcolex<is::cells>(total_energy_density_a);

  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
  auto sie = m.template mdcolex<is::cells>(sie_a);
  auto sound_speed = m.template mdcolex<is::cells>(sound_speed_a);

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      const int old_mat = old_material(i);
      const int new_mat = new_material(i);

      if(old_mat == new_mat) {
        return;
      }

      assert(
        new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_disp.size());

      const double r_old = rho(i);
      const vec<D> u_old = (1.0 / r_old) * mom(i);
      const double e_old =
        (rhoE(i) - 0.5 * r_old * u_old.norm_squared()) / r_old;

      const bool equivalent_eos = eos_pair_equivalent_at_cell<D>(
        eos_disp, old_mat, new_mat, r_old, mom(i), rhoE(i));

      if(equivalent_eos) {
        const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

        velocity(i) = u_old;
        sie(i) = e_old;
        pressure(i) = eos_new.pRhoSie(r_old, e_old);
        sound_speed(i) = eos_new.cRhoSie(r_old, e_old);

        return;
      }

      const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
      const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

      const double p_copy = eos_old.pRhoSie(r_old, e_old);
      const vec<D> u_copy = u_old;

      const auto & q = (new_mat == 0) ? q_ext_plus(i) : q_ext_minus(i);
      const double r_new = q.rho();

      assert(std::isfinite(r_new) && r_new > 0.0);
      assert(std::isfinite(p_copy) && p_copy > 0.0);

      const double e_new = hard::tasks::util::find_sie(eos_new, r_new, p_copy);
      const double c_new = eos_new.cRhoSie(r_new, e_new);

      rho(i) = r_new;
      velocity(i) = u_copy;
      pressure(i) = p_copy;
      sie(i) = e_new;
      sound_speed(i) = c_new;

      mom(i) = r_new * u_copy;
      rhoE(i) = r_new * e_new + 0.5 * r_new * u_copy.norm_squared();
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(rho,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;

      const int old_mat = old_material(i, j);
      const int new_mat = new_material(i, j);

      if(old_mat == new_mat) {
        return;
      }

      assert(
        new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_disp.size());

      const double r_old = rho(i, j);
      const vec<D> u_old = (1.0 / r_old) * mom(i, j);
      const double e_old =
        (rhoE(i, j) - 0.5 * r_old * u_old.norm_squared()) / r_old;

      const bool equivalent_eos = eos_pair_equivalent_at_cell<D>(
        eos_disp, old_mat, new_mat, r_old, mom(i, j), rhoE(i, j));

      if(equivalent_eos) {
        const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

        velocity(i, j) = u_old;
        sie(i, j) = e_old;
        pressure(i, j) = eos_new.pRhoSie(r_old, e_old);
        sound_speed(i, j) = eos_new.cRhoSie(r_old, e_old);

        return;
      }

      const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
      const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

      const double p_copy = eos_old.pRhoSie(r_old, e_old);
      const vec<D> u_copy = u_old;

      const auto & q = (new_mat == 0) ? q_ext_plus(i, j) : q_ext_minus(i, j);
      const double r_new = q.rho();

      assert(std::isfinite(r_new) && r_new > 0.0);
      assert(std::isfinite(p_copy) && p_copy > 0.0);

      const double e_new = hard::tasks::util::find_sie(eos_new, r_new, p_copy);
      const double c_new = eos_new.cRhoSie(r_new, e_new);

      rho(i, j) = r_new;
      velocity(i, j) = u_copy;
      pressure(i, j) = p_copy;
      sie(i, j) = e_new;
      sound_speed(i, j) = c_new;

      mom(i, j) = r_new * u_copy;
      rhoE(i, j) = r_new * e_new + 0.5 * r_new * u_copy.norm_squared();
    };
  }
  else {
    auto pol = get_mdiota_policy(rho,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;

      const int old_mat = old_material(i, j, k);
      const int new_mat = new_material(i, j, k);

      if(old_mat == new_mat) {
        return;
      }

      assert(
        new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_disp.size());

      const double r_old = rho(i, j, k);
      const vec<D> u_old = (1.0 / r_old) * mom(i, j, k);
      const double e_old =
        (rhoE(i, j, k) - 0.5 * r_old * u_old.norm_squared()) / r_old;

      const bool equivalent_eos = eos_pair_equivalent_at_cell<D>(
        eos_disp, old_mat, new_mat, r_old, mom(i, j, k), rhoE(i, j, k));

      if(equivalent_eos) {
        const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

        velocity(i, j, k) = u_old;
        sie(i, j, k) = e_old;
        pressure(i, j, k) = eos_new.pRhoSie(r_old, e_old);
        sound_speed(i, j, k) = eos_new.cRhoSie(r_old, e_old);

        return;
      }

      const auto & eos_old = eos_disp[static_cast<std::size_t>(old_mat)];
      const auto & eos_new = eos_disp[static_cast<std::size_t>(new_mat)];

      const double p_copy = eos_old.pRhoSie(r_old, e_old);
      const vec<D> u_copy = u_old;

      const auto & q =
        (new_mat == 0) ? q_ext_plus(i, j, k) : q_ext_minus(i, j, k);
      const double r_new = q.rho();

      assert(std::isfinite(r_new) && r_new > 0.0);
      assert(std::isfinite(p_copy) && p_copy > 0.0);

      const double e_new = hard::tasks::util::find_sie(eos_new, r_new, p_copy);
      const double c_new = eos_new.cRhoSie(r_new, e_new);

      rho(i, j, k) = r_new;
      velocity(i, j, k) = u_copy;
      pressure(i, j, k) = p_copy;
      sie(i, j, k) = e_new;
      sound_speed(i, j, k) = c_new;

      mom(i, j, k) = r_new * u_copy;
      rhoE(i, j, k) = r_new * e_new + 0.5 * r_new * u_copy.norm_squared();
    };
  }
}

// -----------------------------------------------------------------------------
// 1D diagnostic: print cells whose material changed during determine_psi.
// -----------------------------------------------------------------------------
template<std::size_t D>
void
diagnose_material_changes_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> old_material_a,
  typename field<int>::template accessor<ro, ro> new_material_a,
  field<double>::accessor<ro, ro> mass_density_a,
  typename field<vec<D>>::template accessor<ro, ro> momentum_density_a,
  field<double>::accessor<ro, ro> total_energy_density_a,
  const std::vector<eos::eos_wrapper> & eos_mm,
  int max_print) noexcept {

  if constexpr(D == 1) {
    auto old_material = m.template mdcolex<is::cells>(old_material_a);
    auto new_material = m.template mdcolex<is::cells>(new_material_a);
    auto rho = m.template mdcolex<is::cells>(mass_density_a);
    auto mom = m.template mdcolex<is::cells>(momentum_density_a);
    auto rhoE = m.template mdcolex<is::cells>(total_energy_density_a);

    int nchanged = 0;
    int nprinted = 0;

    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const int old_mat = old_material(i);
      const int new_mat = new_material(i);

      if(old_mat == new_mat) {
        continue;
      }

      ++nchanged;

      const double r = rho(i);
      const double u = mom(i)[0] / r;
      const double e = (rhoE(i) - 0.5 * r * u * u) / r;

      double p_old = -1.0;
      double c_old = -1.0;
      double p_new = -1.0;
      double c_new = -1.0;

      if(old_mat >= 0 && static_cast<std::size_t>(old_mat) < eos_mm.size()) {
        const auto & eos_old = eos_mm[static_cast<std::size_t>(old_mat)];
        p_old = eos_old.pRhoSie(r, e);
        c_old = eos_old.cRhoSie(r, e);
      }

      if(new_mat >= 0 && static_cast<std::size_t>(new_mat) < eos_mm.size()) {
        const auto & eos_new = eos_mm[static_cast<std::size_t>(new_mat)];
        p_new = eos_new.pRhoSie(r, e);
        c_new = eos_new.cRhoSie(r, e);
      }

      if(nprinted < max_print) {
        std::cout << "[material changed] i=" << i << " old=" << old_mat
                  << " new=" << new_mat << " rho=" << r << " u=" << u
                  << " rhoE=" << rhoE(i) << " e_from_cons=" << e
                  << " p_oldEOS=" << p_old << " c_oldEOS=" << c_old
                  << " p_newEOS=" << p_new << " c_newEOS=" << c_new << "\n";

        ++nprinted;
      }
    }

    std::cout << "[material changed total] " << nchanged << "\n";
  }
}

struct diagnose_psi_range_1d_task : annotation_task_region {
  inline static const std::string name{"diagnose_psi_range_1d"};
};

template<std::size_t D>
void
diagnose_psi_range_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, ro> psi_a,
  typename field<vec<D>>::template accessor<ro, ro> velocity_a,
  const char * label) noexcept {

  if constexpr(D == 1) {
    auto psi = m.template mdcolex<is::cells>(psi_a);
    auto vel = m.template mdcolex<is::cells>(velocity_a);

    double min_psi = std::numeric_limits<double>::max();
    double max_psi = -std::numeric_limits<double>::max();
    double max_abs_psi = 0.0;
    double min_u = std::numeric_limits<double>::max();
    double max_u = -std::numeric_limits<double>::max();

    int bad = 0;

    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const double p = psi(i);
      const double u = vel(i)[0];

      min_psi = std::min(min_psi, p);
      max_psi = std::max(max_psi, p);
      max_abs_psi = std::max(max_abs_psi, std::abs(p));
      min_u = std::min(min_u, u);
      max_u = std::max(max_u, u);

      if(!std::isfinite(p) || std::abs(p) > 10.0) {
        ++bad;
        if(bad <= 10) {
          std::cout << "[bad psi] " << label << " i=" << i << " psi=" << p
                    << " u=" << u << "\n";
        }
      }
    }

    if(bad > 0) {
      std::cout << "[psi range] " << label << " min=" << min_psi
                << " max=" << max_psi << " maxabs=" << max_abs_psi
                << " u_min=" << min_u << " u_max=" << max_u
                << " bad_count=" << bad << "\n";
    }
  }
}

template<std::size_t D>
void
diagnose_material_cuts_1d(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> material_a,
  field<double>::accessor<ro, ro> psi_a) noexcept {

  if constexpr(D == 1) {
    auto material = m.template mdcolex<is::cells>(material_a);
    auto psi = m.template mdcolex<is::cells>(psi_a);

    int ncuts = 0;

    for(auto i : m.template cells<ax::x, dm::corrector>()) {
      if(material(i - 1) != material(i)) {
        ++ncuts;

        const double xL = 0.5 * (m.template left<ax::x>(i - 1) +
                                  m.template right<ax::x>(i - 1));
        const double xR =
          0.5 * (m.template left<ax::x>(i) + m.template right<ax::x>(i));

        std::cout << "[material cut] i=" << i << " xL=" << xL << " xR=" << xR
                  << " mL=" << material(i - 1) << " mR=" << material(i)
                  << " psiL=" << psi(i - 1) << " psiR=" << psi(i) << "\n";
      }
    }

    std::cout << "[material cuts total] " << ncuts << "\n";
  }
}

} // namespace hard::tasks::level_set
#endif
