#define IN_TEST

#include "flecsi/execution.hh"
#include "flecsi/util/unit.hh"

#include "spec/labels.hh"

// NOTE: These #includes need to keep this relative order
// clang-format off
#include "tests_level_set_utils.hh"
#include "../../../modules/gfm/tasks/ls_core.hh"
// clang-format on

using namespace flecsi::exec;

using namespace hard::tasks::level_set;
using namespace tests::level_set_utils;
using flecsi::execute;
using flecsi::scheduler;

template<std::size_t D>
typename mesh<D>::periodic_axes
init_boundaries() noexcept {
  typename mesh<D>::periodic_axes p;
  p[ax::x] = false;
  if constexpr(D >= 2) {
    p[ax::y] = false;
  }
  if constexpr(D == 3) {
    p[ax::z] = false;
  } // if
  return p;
}

template<std::size_t D>
void
init_mesh(scheduler & s, state<D> & st) {

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

  std::size_t levels = 5;

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
init_dt(flecsi::exec::accelerator s,
  typename single<double>::template accessor<wo> t,
  typename single<double>::template accessor<wo> dt) noexcept {
  *t = 0;
  *dt = 0.1; // TODO: Make this work with CFL = 0.4 (CALCULATE)
}

template<std::size_t D>
void
init_fields(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<wo, wo> psi_a,
  typename field<int>::template accessor<wo, wo> material_a) noexcept {

  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto mat = m.template mdcolex<is::cells>(material_a);

  // If your domain is [0,1]^D, these are normalized.
  // If not, scale x0,y0,z0,R accordingly.
  const double x0 = 0.5;
  const double y0 = 0.75;
  const double z0 = 0.5;
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
advance_t(flecsi::exec::accelerator s,
  typename single<double>::template accessor<rw> t,
  typename single<double>::template accessor<ro> dt) noexcept {
  *t += *dt;
}; // advance_t

template<std::size_t D>
double
compare_material(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, na> material_a,
  typename field<int>::template accessor<ro, na> material_comp_a) noexcept {

  auto mat = m.template mdcolex<is::cells>(material_a);
  auto matc = m.template mdcolex<is::cells>(material_comp_a);

  if constexpr(D == 1) {
    auto dx = m.template delta<ax::x>();
    return std::sqrt(s.executor().reduceall(
      i, l2, (m.template cells<ax::x, dm::quantities>()), fold::sum, double) {
      l2(std::pow((mat(i) - matc(i)) * dx, 2));
    });
  }
  else if constexpr(D == 2) {
    auto dx = m.template delta<ax::x>();
    auto dy = m.template delta<ax::y>();
    auto pol = get_mdiota_policy(mat,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return std::sqrt(s.executor().reduceall(ji, l2, pol, fold::sum, double) {
      auto [j, i] = ji;
      l2(std::pow((mat(i, j) - matc(i, j)) * dx * dy, 2));
    });
  }
  else { // D == 3
    auto dx = m.template delta<ax::x>();
    auto dy = m.template delta<ax::y>();
    auto dz = m.template delta<ax::z>();
    auto pol = get_mdiota_policy(mat,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return std::sqrt(s.executor().reduceall(kji, l2, pol, fold::sum, double) {
      auto [k, j, i] = kji;
      l2(std::pow((mat(i, j, k) - matc(i, j, k)) * dx * dy * dz, 2));
    });
  }
} // compare_material

inline void
print_time_dt(flecsi::exec::accelerator,
  typename single<double>::template accessor<ro> t,
  typename single<double>::template accessor<ro> dt) noexcept {
  std::cout << "[debug] t=" << *t << " dt=" << *dt << std::endl;
}

template<std::size_t D>
double
do_test(flecsi::scheduler & s, state<D> & st) {

  s.execute<copy_scal_field<dm::quantities, D, int>>(
    flecsi::exec::on, *st.m, st.material(*st.m), st.initial_material(*st.m));

  // TOOD : remove below
  // ---- 1) Print initial time/dt (must be via a task) ----
  s.execute<print_time_dt>(flecsi::exec::on, st.t(*st.gt), st.dt(*st.gt));

  // ---- 2) Choose end time / step count ----
  const double t_end = 10; // <-- change this
  const double dt_host =
    0.1; // must match init_dt unless you add a "read dt" task
  const std::size_t nsteps =
    static_cast<std::size_t>(std::lround(t_end / dt_host));

  std::cout << "[debug] nsteps=" << nsteps << " (t_end=" << t_end
            << ", dt=" << dt_host << ")\n";

  // This test does a time loop (t: 0 -> 2T). At the end, the system should
  // have returned to its original state
  for(std::size_t n{0}; n < nsteps; n++) {
    ssprk3_step(s, st); // uses rhs + frozen_material internally
    s.execute<advance_t>(flecsi::exec::on, st.t(*st.gt), st.dt(*st.gt));
  }
  // ---- 4) Print final time/dt ----
  s.execute<print_time_dt>(flecsi::exec::on, st.t(*st.gt), st.dt(*st.gt));

  // TODO: remove belows
  s.execute<print_scal_field<dm::quantities, D, int>>(
    flecsi::exec::on, *st.m, st.initial_material(*st.m));
  s.execute<print_scal_field<dm::quantities, D, int>>(
    flecsi::exec::on, *st.m, st.material(*st.m));

  return s
    .reduce<compare_material<D>, flecsi::exec::fold::sum>(
      flecsi::exec::on, *st.m, st.material(*st.m), st.initial_material(*st.m))
    .get();
}

int
test_driver(scheduler & s) {
  UNIT() {
    // TODO: Uncomment
    // // 1D
    // state<1> d1;
    // init_mesh<1>(s, d1);
    // s.execute<init_fields<1>>(
    //   flecsi::exec::on, *d1.m, d1.psi(*d1.m), d1.label(*d1.m));
    // s.execute<init_dt>(flecsi::exec::on, d1.t(*d1.gt), d1.dt(*d1.gt));
    // EXPECT_LT(do_test<1>(s, d1), 1e-10);

    // 2D
    state<2> d2;
    init_mesh<2>(s, d2);
    s.execute<init_fields<2>>(
      flecsi::exec::on, *d2.m, d2.psi(*d2.m), d2.material(*d2.m));
    s.execute<init_dt>(flecsi::exec::on, d2.t(*d2.gt), d2.dt(*d2.gt));
    EXPECT_LT(do_test<2>(s, d2), 1e-10);

    // TODO: Uncomment
    // // 3D
    // state<3> d3;
    // init_mesh<3>(s, d3);
    // s.execute<init_fields<3>>(
    //   flecsi::exec::on, *d3.m, d3.psi(*d3.m), d3.label(*d3.m));
    // s.execute<init_dt>(flecsi::exec::on, d3.t(*d3.gt), d3.dt(*d3.gt));
    // EXPECT_LT(do_test<3>(s, d3), 1e-10);
  };
} // kernel_driver

flecsi::util::unit::driver<test_driver> driver;
