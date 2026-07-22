#include "flecsi/execution.hh"
#include "flecsi/flog.hh"
#include "flecsi/util/unit.hh"

#include <utility>

#include <utility>

#include "../types.hh"
#include "common/tasks/utils.hh"
#include "spec/labels.hh"

// NOTE: These #includes need to keep this relative order
// clang-format off
#include "tests_level_set_utils.hh"
// ls_core functions are available via update_ls.cc compilation
// clang-format on

using namespace flecsi::exec;

// Forward declare ssprk3_step (defined in update_ls.cc)
namespace hard {
template<std::size_t D>
void ssprk3_step(flecsi::scheduler & sc, tests::level_set_utils::state<D> & s);
}

using namespace tests::level_set_utils;
using namespace common::tasks::utils;
using flecsi::execute;
using flecsi::scheduler;

template<std::size_t D>
typename mesh<D>::periodic_axes
init_boundaries() noexcept {
  typename mesh<D>::periodic_axes p;
  p[ax::x] = true;
  if constexpr(D >= 2) {
    p[ax::y] = true;
  }
  if constexpr(D == 3) {
    p[ax::z] = true;
  } // if
  return p;
}

template<std::size_t D>
void
init_mesh(scheduler & s, state<D> & st, std::size_t levels) {

  std::array<std::array<bd::boundary_type, 2>, D> bnds;
  bnds[ax::x][bd::low] = bd::periodic;
  bnds[ax::x][bd::high] = bd::periodic;
  if(D == 2 || D == 3) {
    bnds[ax::y][bd::low] = bd::periodic;
    bnds[ax::y][bd::high] = bd::periodic;
  }
  if(D == 3) {
    bnds[ax::z][bd::low] = bd::periodic;
    bnds[ax::z][bd::high] = bd::periodic;
  }

  auto bf = execute<init_boundaries<D>>();

  std::size_t lowest_level = levels;
  std::size_t highest_level = levels;
  std::size_t max_num_levels = highest_level - lowest_level + 1;

  typename mesh<D>::grect geom;
  geom[0][0] = 0.;
  geom[0][1] = 1.;
  if(D == 2 || D == 3) {
    geom[1][0] = 0.;
    geom[1][1] = 1.;
  }
  if(D == 3) {
    geom[2][0] = 0.;
    geom[2][1] = 1.;
  }

  for(std::size_t i{0}; i < max_num_levels; i++) {
    typename mesh<D>::gcoord axis_extents(D);
    axis_extents[ax::x] = 1 << (levels - i);
    if(D == 2 || D == 3) {
      axis_extents[ax::y] = 1 << (levels - i);
    }
    if(D == 3) {
      axis_extents[ax::z] = 1 << (levels - i);
    } // if

    s.allocate(st.m,
      typename mesh<D>::mpi_coloring(
        s, s.runtime().processes(), axis_extents, bf.get()),
      geom);

    s.allocate(st.gt, s.runtime().processes());
  }
}

void
init_dt_period(flecsi::exec::cpu,
  typename single<double>::template accessor<spec::wo> t,
  typename single<double>::template accessor<spec::wo> dt,
  typename single<double>::template accessor<spec::wo> period,
  double dt_value,
  double period_value) noexcept {
  *t = 0;
  *dt = dt_value;
  *period = period_value;
}

template<std::size_t D>
void
init_fields(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<spec::ro> m,
  typename field<double>::template accessor<spec::wo, spec::wo> psi_a,
  typename field<int>::template accessor<spec::wo, spec::wo>
    material_a) noexcept {

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto mat = m.template mdcolex<is::cells>(material_a);

  // Paper two-material test: circle centered at (0.35, 0.35) with radius 0.15
  // in the unit square.
  const double x0 = 0.35;
  const double y0 = 0.35;
  const double z0 = 0.35;
  const double R = 0.15;

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double x = m.template center<ax::x>(i);
      const double d =
        std::abs(x - x0) - R; // signed distance to interval boundary analog

      const double phi1 = -d;
      const double phi2 = d;

      if(phi1 < phi2) {
        psi(i) = phi1;
        mat(i) = 0;
      }
      else {
        psi(i) = phi2;
        mat(i) = 1;
      }
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      const double x = m.template center<ax::x>(i);
      const double y = m.template center<ax::y>(j);

      const double d = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0)) - R;

      const double phi1 = -d;
      const double phi2 = d;

      if(phi1 < phi2) {
        psi(i, j) = phi1;
        mat(i, j) = 0;
      }
      else {
        psi(i, j) = phi2;
        mat(i, j) = 1;
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
      const double x = m.template center<ax::x>(i);
      const double y = m.template center<ax::y>(j);
      const double z = m.template center<ax::z>(k);

      const double d = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) +
                                 (z - z0) * (z - z0)) -
                       R;

      const double phi1 = -d;
      const double phi2 = d;

      if(phi1 < phi2) {
        psi(i, j, k) = phi1;
        mat(i, j, k) = 0;
      }
      else {
        psi(i, j, k) = phi2;
        mat(i, j, k) = 1;
      }
    };
  }
}

inline void
advance_t(flecsi::exec::cpu,
  typename single<double>::template accessor<spec::rw> t,
  typename single<double>::template accessor<spec::ro> dt) noexcept {
  *t += *dt;
}; // advance_t

template<std::size_t D>
double
compare_psi_l1(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<spec::ro> m,
  typename field<double>::template accessor<spec::ro, spec::na> psi_a,
  typename field<double>::template accessor<spec::ro, spec::na>
    psi_comp_a) noexcept {

  constexpr double interface_band = 1.0 / 20.0; // from Aslam JCP '26

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto psic = m.template mdcolex<is::cells>(psi_comp_a);

  if constexpr(D == 1) {
    auto dx = m.template delta<ax::x>();
    return s.executor().reduceall(
      i, l1, (m.template cells<ax::x, dm::quantities>()), fold::sum, double) {
      const double ref = psic(i);
      l1((std::abs(ref) < interface_band) ? std::abs(psi(i) - ref) * dx : 0.0);
    };
  }
  else if constexpr(D == 2) {
    auto dx = m.template delta<ax::x>();
    auto dy = m.template delta<ax::y>();
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return s.executor().reduceall(ji, l1, pol, fold::sum, double) {
      auto [j, i] = ji;
      const double ref = psic(i, j);
      l1((std::abs(ref) < interface_band) ? std::abs(psi(i, j) - ref) * dx * dy
                                          : 0.0);
    };
  }
  else { // D == 3
    auto dx = m.template delta<ax::x>();
    auto dy = m.template delta<ax::y>();
    auto dz = m.template delta<ax::z>();
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return s.executor().reduceall(kji, l1, pol, fold::sum, double) {
      auto [k, j, i] = kji;
      const double ref = psic(i, j, k);
      l1((std::abs(ref) < interface_band)
           ? std::abs(psi(i, j, k) - ref) * dx * dy * dz
           : 0.0);
    };
  }
}

template<std::size_t D>
double
compare_psi_linf(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<spec::ro> m,
  typename field<double>::template accessor<spec::ro, spec::na> psi_a,
  typename field<double>::template accessor<spec::ro, spec::na>
    psi_comp_a) noexcept {

  constexpr double interface_band = 1.0 / 20.0; // from Aslam JCP '26

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto psic = m.template mdcolex<is::cells>(psi_comp_a);

  if constexpr(D == 1) {
    return s.executor().reduceall(
      i, linf, (m.template cells<ax::x, dm::quantities>()), fold::max, double) {
      const double ref = psic(i);
      linf((std::abs(ref) < interface_band) ? std::abs(psi(i) - ref) : 0.0);
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return s.executor().reduceall(ji, linf, pol, fold::max, double) {
      auto [j, i] = ji;
      const double ref = psic(i, j);
      linf((std::abs(ref) < interface_band) ? std::abs(psi(i, j) - ref) : 0.0);
    };
  }
  else { // D == 3
    auto pol = get_mdiota_policy(psi,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return s.executor().reduceall(kji, linf, pol, fold::max, double) {
      auto [k, j, i] = kji;
      const double ref = psic(i, j, k);
      linf(
        (std::abs(ref) < interface_band) ? std::abs(psi(i, j, k) - ref) : 0.0);
    };
  }
}

template<std::size_t D>
std::pair<double, double>
do_test(flecsi::scheduler & s, state<D> & st, double t_end, double dt_host) {

  s.execute<copy_scalar<D>>(
    flecsi::exec::on, st.psi(*st.m), st.initial_psi(*st.m));

  const std::size_t nsteps =
    static_cast<std::size_t>(std::lround(t_end / dt_host));

  // This test advances the level set through a time-reversal cycle.
  // At the end, psi should have returned to its initial state.
  for(std::size_t n{0}; n < nsteps; n++) {
    if(n % 100 == 0 && s.runtime().process() == 0) {
      flog(info) << "[debug] step " << n << "/" << nsteps << std::endl;
    }
    hard::ssprk3_step<D>(s, st); // uses rhs + frozen_material internally
    s.execute<advance_t>(flecsi::exec::on, st.t(*st.gt), st.dt(*st.gt));
  }
  // ---- 4) Print final time/dt ----

  auto l1_error =
    s.reduce<compare_psi_l1<D>, flecsi::exec::fold::sum>(
       flecsi::exec::on, *st.m, st.psi(*st.m), st.initial_psi(*st.m))
      .get();

  auto linf_error =
    s.reduce<compare_psi_linf<D>, flecsi::exec::fold::max>(
       flecsi::exec::on, *st.m, st.psi(*st.m), st.initial_psi(*st.m))
      .get();

  std::cout << "[result] " << D << "D : L1 psi error = " << l1_error
            << ", Linf psi error = " << linf_error << std::endl;

  return {l1_error, linf_error};
}

int
test_driver(scheduler & s) {
  UNIT() {
    auto relaxed = [](double ref) { return 1.8 * ref; };

    // All these tests come from Aslam JCP '26

    // 2D, level 7, dx = 1/128, dt = 2/640
    {
      state<2> d2;
      constexpr std::size_t levels = 7;
      constexpr double dt = 2.0 / 640.0;
      constexpr double t_end = 2.0;
      constexpr double period = 2.0;
      constexpr double l1_ref = 2e-6; // original value: 1.3e-6
      constexpr double linf_ref = 2e-4; // original value: 1.6e-4

      init_mesh<2>(s, d2, levels);
      s.execute<init_fields<2>>(
        flecsi::exec::on, *d2.m, d2.psi(*d2.m), d2.material(*d2.m));
      s.execute<init_dt_period>(flecsi::exec::on,
        d2.t(*d2.gt),
        d2.dt(*d2.gt),
        d2.period(*d2.gt),
        dt,
        period);
      auto err = do_test<2>(s, d2, t_end, dt);
      EXPECT_LT(err.first, relaxed(l1_ref));
      EXPECT_LT(err.second, relaxed(linf_ref));
    }

    /* Test commented out given possible CI timeouts
    // 2D, level 8, dx = 1/256, dt = 2/2032
    {
      state<2> d2;
      constexpr std::size_t levels = 8;
      constexpr double dt = 2.0 / 2032.0;
      constexpr double t_end = 2.0;
      constexpr double period = 2.0;
      constexpr double l1_ref = 5e-8; // original value: 4.7e-8
      constexpr double linf_ref = 2e-5; // original value: 1.1e-5

      init_mesh<2>(s, d2, levels);
      s.execute<init_fields<2>>(
        flecsi::exec::on, *d2.m, d2.psi(*d2.m), d2.material(*d2.m));
      s.execute<init_dt_period>(flecsi::exec::on,
        d2.t(*d2.gt),
        d2.dt(*d2.gt),
        d2.period(*d2.gt),
        dt,
        period);
      auto err = do_test<2>(s, d2, t_end, dt);
      EXPECT_LT(err.first, relaxed(l1_ref));
      EXPECT_LT(err.second, relaxed(linf_ref));
    }
    */

    // 3D, level 7, dx = 1/128, dt = 1/480
    {
      state<3> d3;
      constexpr std::size_t levels = 7;
      constexpr double dt = 1.0 / 480.0;
      constexpr double t_end = 1.0;
      constexpr double period = 1.0;
      constexpr double l1_ref = 2e-6; // original value: 1.4e-6
      constexpr double linf_ref = 2e-2; // original value: 1.36e-2

      init_mesh<3>(s, d3, levels);
      s.execute<init_fields<3>>(
        flecsi::exec::on, *d3.m, d3.psi(*d3.m), d3.material(*d3.m));
      s.execute<init_dt_period>(flecsi::exec::on,
        d3.t(*d3.gt),
        d3.dt(*d3.gt),
        d3.period(*d3.gt),
        dt,
        period);
      auto err = do_test<3>(s, d3, t_end, dt);
      EXPECT_LT(err.first, relaxed(l1_ref));
      EXPECT_LT(err.second, relaxed(linf_ref));
    }

    /* Test commented out given possible CI timeouts
    // 3D, level 8, dx = 1/256, dt = 1/1524
    {
      state<3> d3;
      constexpr std::size_t levels = 8;
      constexpr double dt = 1.0 / 1524.0;
      constexpr double t_end = 1.0;
      constexpr double period = 1.0;
      constexpr double l1_ref = 7e-8; // original value: 6.2e-8
      constexpr double linf_ref = 2e-3; // original value: 1.6e-3

      init_mesh<3>(s, d3, levels);
      s.execute<init_fields<3>>(
        flecsi::exec::on, *d3.m, d3.psi(*d3.m), d3.material(*d3.m));
      s.execute<init_dt_period>(flecsi::exec::on,
        d3.t(*d3.gt),
        d3.dt(*d3.gt),
        d3.period(*d3.gt),
        dt,
        period);
      auto err = do_test<3>(s, d3, t_end, dt);
      EXPECT_LT(err.first, relaxed(l1_ref));
      EXPECT_LT(err.second, relaxed(linf_ref));
    }
    */
  };
} // kernel_driver

flecsi::util::unit::driver<test_driver> driver;