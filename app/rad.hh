#ifndef HARD_RAD_HH
#define HARD_RAD_HH

#include "state.hh"
#include "tasks/rad.hh"

#include <cstddef>

namespace hard::rad {

template<std::size_t D>
void
vcycle(control_policy<state, D> & cp, std::size_t index) {
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
      s.mgr.Esf.flip();
      // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
      // parallelized
      sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
        mf,
        s.mgr.Ew(mf),
        s.mgr.Esf(mf),
        s.mgr.Esf(mf, 1),
        s.mgr.Ef(mf),
        0.8);
    } // for
  }
  else {

    // FIXME: Remove when finished with debugging
    // flog(warn) << "Cycle level(index): " << level << "(" << index << ")"
    //            << std::endl;

    auto & mc = *s.mh[index + 1];

    // Pre Smoothing
    for(std::size_t i{0}; i < s.mgr.mg_pre; ++i) {
      s.mgr.Esf.flip();
      sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
        mf,
        s.mgr.Ew(mf),
        s.mgr.Esf(mf),
        s.mgr.Esf(mf, 1),
        s.mgr.Ef(mf),
        0.8);
    } // for

    // Recursive solve
    sc.execute<task::rad::residual<D>>(flecsi::exec::on,
      mf,
      s.mgr.Ew(mf),
      s.mgr.Esf(mf),
      s.mgr.Ef(mf),
      s.mgr.Resf(mf));

    sc.execute<task::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.mgr.Resf(mf), s.mgr.Ef(mc));

    // Initialize the solution fields for the coarser level
    sc.execute<task::rad::const_init<D>>(
      flecsi::exec::on, mc, s.mgr.Esf(mc), 0.0);
    sc.execute<task::rad::const_init<D>>(
      flecsi::exec::on, mc, s.mgr.Esf(mc, 1), 0.0);

    vcycle<D>(cp, index + 1);

    sc.execute<task::rad::nlinear_interpolation<D>>(
      flecsi::exec::on, mc, mf, s.mgr.Esf(mc), s.mgr.Errf(mf));

    sc.execute<task::rad::correction<D>>(
      flecsi::exec::on, mf, s.mgr.Esf(mf), s.mgr.Errf(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < s.mgr.mg_post; ++i) {
      s.mgr.Esf.flip();
      sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
        mf,
        s.mgr.Ew(mf),
        s.mgr.Esf(mf),
        s.mgr.Esf(mf, 1),
        s.mgr.Ef(mf),
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
    sc.execute<task::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.mgr.Ef(mf), s.mgr.Ef(mc));
    sc.execute<task::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.mgr.Esf(mf), s.mgr.Esf(mc));

    // Set the diffusion coefficient and the stencil (TODO)
    sc.execute<task::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.mgr.Df_x(mf), s.mgr.Df_x(mc));

    sc.execute<task::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.mgr.Df_y(mf), s.mgr.Df_y(mc));

    sc.execute<task::rad::full_weighting<D>>(
      flecsi::exec::on, mf, mc, s.mgr.Df_z(mf), s.mgr.Df_z(mc));

    sc.execute<task::rad::stencil_init<D>>(flecsi::exec::on,
      mc,
      s.mgr.Df_x(mc),
      s.mgr.Df_y(mc),
      s.mgr.Df_z(mc),
      s.mgr.Ew(mc),
      s.dt(*s.gt));

    // Now call solve for one level deeper
    fmg<D>(cp, index + 1);

    // Interpolate solution back up (RHS does not change)
    sc.execute<task::rad::nlinear_interpolation<D>>(
      flecsi::exec::on, mc, mf, s.mgr.Esf(mc), s.mgr.Esf(mf));

    // Do a V-Cycle
    for(std::size_t i{0}; i < s.mgr.mg_cycles; ++i) {
      vcycle<D>(cp, index);
    } // for
  } // if
} // fmg

} // namespace hard::rad

#endif // HARD_RAD_HH
