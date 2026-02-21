#ifndef HARD_MODULES_RAD_LINSOLVE_HH
#define HARD_MODULES_RAD_LINSOLVE_HH

#include "flecsi/execution.hh"
#include "flecsi/flog.hh"

#include "flecsolvers.hh"
#include "tasks/rad.hh"

namespace hard {

template<std::size_t D>
void
vcycle(control_policy<state, D> & cp, std::size_t index) {

  using namespace rad;

  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();
  auto & mf = *s.mh[index];

  // Find current level
  std::size_t level{s.highest_level - index};

  if(level == s.lowest_level) {

    // FIXME: Remove when finished with debugging
    // flog(warn) << "Direct solve level(index): " << level << "(" << index <<
    // ")"
    //            << std::endl;

    // Direct solve for a single interior point
    for(std::size_t i{0}; i < 100; i++) {
      s.rad.mgr.Esf.flip();
      // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
      // parallelized
      sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
        mf,
        s.rad.mgr.Ew(mf),
        s.rad.mgr.Esf(mf),
        s.rad.mgr.Esf(mf, 1),
        s.rad.mgr.Ef(mf),
        0.8);
    } // for
  }
  else {

    // FIXME: Remove when finished with debugging
    // flog(warn) << "Cycle level(index): " << level << "(" << index << ")"
    //            << std::endl;

    auto & mc = *s.mh[index + 1];

    // Pre Smoothing
    for(std::size_t i{0}; i < s.rad.mgr.mg_pre; ++i) {
      s.rad.mgr.Esf.flip();
      sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
        mf,
        s.rad.mgr.Ew(mf),
        s.rad.mgr.Esf(mf),
        s.rad.mgr.Esf(mf, 1),
        s.rad.mgr.Ef(mf),
        0.8);
    } // for

    // Recursive solve
    sc.execute<tasks::rad::residual<D>>(flecsi::exec::on,
      mf,
      s.rad.mgr.Ew(mf),
      s.rad.mgr.Esf(mf),
      s.rad.mgr.Ef(mf),
      s.rad.mgr.Resf(mf));

    sc.execute<tasks::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.rad.mgr.Resf(mf), s.rad.mgr.Ef(mc));

    // Initialize the solution fields for the coarser level
    sc.execute<tasks::rad::const_init<D>>(
      flecsi::exec::on, mc, s.rad.mgr.Esf(mc), 0.0);
    sc.execute<tasks::rad::const_init<D>>(
      flecsi::exec::on, mc, s.rad.mgr.Esf(mc, 1), 0.0);

    vcycle<D>(cp, index + 1);

    sc.execute<tasks::rad::nlinear_interpolation<D>>(
      flecsi::exec::on, mc, mf, s.rad.mgr.Esf(mc), s.rad.mgr.Errf(mf));

    sc.execute<tasks::rad::correction<D>>(
      flecsi::exec::on, mf, s.rad.mgr.Esf(mf), s.rad.mgr.Errf(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < s.rad.mgr.mg_post; ++i) {
      s.rad.mgr.Esf.flip();
      sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
        mf,
        s.rad.mgr.Ew(mf),
        s.rad.mgr.Esf(mf),
        s.rad.mgr.Esf(mf, 1),
        s.rad.mgr.Ef(mf),
        0.8);
    } // for
  } // if
} // vcycle

template<std::size_t D>
void
fmg(control_policy<state, D> & cp, std::size_t index = 0) {
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();
  auto & mf = *s.mh[index];

  // The scheme requires:
  // 1) Go to a coarser grid, adapt all and repeat this step
  // 2) If in the deeper level, direct solve or do a V-Cycle for a number of

  // iterations 3) Come back up, interpolate, and do a V-Cycle

  // Find current level
  std::size_t level{s.highest_level - index};

  // Deepest level
  if(level == s.lowest_level) {

    // FIXME: Remove when finished with debugging
    // flog(warn) << "Deepest level(index):" << level << "(" << index << ")"
    //            << std::endl;

    // If in the deepest level, the V-Cycle is already doing a direct solve
    vcycle<D>(cp, index);
  }
  else {

    // FIXME: Remove when finished with debugging
    // flog(warn) << "cycle level(index): " << level << "(" << index << ")"
    //            << std::endl;
    auto & mc = *s.mh[index + 1];

    // Set the RHS and solution field
    sc.execute<tasks::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.rad.mgr.Ef(mf), s.rad.mgr.Ef(mc));
    sc.execute<tasks::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.rad.mgr.Esf(mf), s.rad.mgr.Esf(mc));

    // Set the diffusion coefficient and the stencil (TODO)
    sc.execute<tasks::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.rad.mgr.Df_x(mf), s.rad.mgr.Df_x(mc));

    sc.execute<tasks::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.rad.mgr.Df_y(mf), s.rad.mgr.Df_y(mc));

    sc.execute<tasks::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.rad.mgr.Df_z(mf), s.rad.mgr.Df_z(mc));

    sc.execute<tasks::rad::stencil_init<D>>(flecsi::exec::on,
      mc,
      s.rad.mgr.Df_x(mc),
      s.rad.mgr.Df_y(mc),
      s.rad.mgr.Df_z(mc),
      s.rad.mgr.Ew(mc),
      s.dt(*s.gt));

    // Now call solve for one level deeper
    fmg<D>(cp, index + 1);

    // Interpolate solution back up (RHS does not change)
    sc.execute<tasks::rad::nlinear_interpolation<D>>(
      flecsi::exec::on, mc, mf, s.rad.mgr.Esf(mc), s.rad.mgr.Esf(mf));

    // Do a V-Cycle
    for(std::size_t i{0}; i < s.rad.mgr.mg_cycles; ++i) {
      vcycle<D>(cp, index);
    } // for
  } // if
} // fmg

template<std::size_t D>
auto
make_solver(control_policy<state, D> & cp) {

  auto & s = cp.state();
  using namespace flecsolve;
  using namespace rad;

  // Solver parameters
  double temp{0.0};

  // Only pass cp instead?
  solver_parameters<D> params{
    std::ref(cp.state()), std::ref(cp.scheduler()), temp};
  op::core<operator_t<D>> so(params);
  auto op_handle = op::ref(so);

  precond_parameters<D> pparams{std::ref(s),
    std::ref(cp.scheduler()),
    0,
    s.rad.mgr.jacobi_iterations,
    s.rad.mgr.nr_vcycles};
  op::core<v_cycle<D>> po(pparams);
  auto prec_handle = op::ref(po);

  auto f = flecsolve::vec::make(s.rad.mgr.Ef(*s.m));

  auto slv = flecsolve::bicgstab::solver(
    s.rad.mgr.solver_settings, flecsolve::bicgstab::make_work(f))(
    op_handle, prec_handle, [&](auto &, double) { return false; });
  return slv;
}

template<std::size_t D>
void
linsolve(control_policy<state, D> & cp) {

  flecsi::scheduler & sc = cp.scheduler();
  auto & s = cp.state();
  auto & mf = *s.m;

  if(s.rad.mgr.full_multigrid) {
    fmg<D>(cp);

    sc.execute<tasks::rad::residual<D>>(flecsi::exec::on,
      mf,
      s.rad.mgr.Ew(mf),
      s.rad.mgr.Esf(mf),
      s.rad.mgr.Ef(mf),
      s.rad.mgr.Resf(mf));
    auto r = flecsolve::vec::make(s.rad.mgr.Resf(mf));
    flog(info) << "final res norm radiation: " << r.l2norm().get() << std::endl;
  }
  else {
    // flecsolve vectors
    auto f = flecsolve::vec::make(s.rad.mgr.Ef(mf));
    auto u = flecsolve::vec::make(s.rad.mgr.Uf(mf));

    /* Sets up a new linear solver at each iteration (this should not be
     * done).*/
    auto slv = make_solver<D>(cp);
    auto slv_info = slv(f, u);
    auto iters = slv_info.iters;
    auto res_norm_final = slv_info.res_norm_final;

    flog(info) << "final res norm radiation (flecsolve): " << res_norm_final
               << " iter: " << iters << std::endl;

    sc.execute<tasks::rad::residual<D>>(flecsi::exec::on,
      mf,
      s.rad.mgr.Ew(mf),
      s.rad.mgr.Uf(mf),
      s.rad.mgr.Ef(mf),
      s.rad.mgr.Resf(mf));
    auto r = flecsolve::vec::make(s.rad.mgr.Resf(mf));
    flog(info) << "final res norm radiation: " << r.l2norm().get() << std::endl;
  }
} // linsolve
} // namespace hard

#endif
