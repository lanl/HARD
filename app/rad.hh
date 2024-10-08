#ifndef HARD_RAD_HH
#define HARD_RAD_HH

#include "state.hh"
#include "tasks/rad.hh"

#include <cstddef>

namespace hard::rad {

template<std::size_t D>
void
vcycle(state<D> & s, std::size_t index) {
  auto & mf = *s.mh[index];

  // Find current level
  std::size_t level{s.highest_level - index};

  if(level == s.lowest_level) {
#ifndef HARD_BENCHMARK_MODE

    flog(warn) << "Direct solve level(index): " << level << "(" << index << ")"
               << std::endl;

#endif
    // Direct solve for a single interior point
    for(std::size_t i{0}; i < 100; i++) {
      s.Esf.flip();
      // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
      // parallelized
      flecsi::execute<task::rad::damped_jacobi<D>, flecsi::default_accelerator>(
        mf, s.Ew(mf), s.Esf(mf), s.Esf(mf, 1), s.Ef(mf), 0.8);
      // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf));
      // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf, 1));
    } // for
  }
  else {
#ifndef HARD_BENCHMARK_MODE

    flog(warn) << "Cycle level(index): " << level << "(" << index << ")"
               << std::endl;

#endif
    auto & mc = *s.mh[index + 1];

    // Pre Smoothing
    for(std::size_t i{0}; i < s.mg_pre; ++i) {
      s.Esf.flip();
      flecsi::execute<task::rad::damped_jacobi<D>, flecsi::default_accelerator>(
        mf, s.Ew(mf), s.Esf(mf), s.Esf(mf, 1), s.Ef(mf), 0.8);
      // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf));
      // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf, 1));
    } // for

    // Recursive solve
    flecsi::execute<task::rad::residual<D>, flecsi::default_accelerator>(
      mf, s.Ew(mf), s.Esf(mf), s.Ef(mf), s.Resf(mf));
    // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Resf(mf));

    flecsi::execute<task::rad::full_weighting<D>, flecsi::default_accelerator>(
      mf, mc, s.Resf(mf), s.Ef(mc));
    // flecsi::execute<task::rad::apply_BC<D>>(mc, s.Ef(mc));

    // Initialize the solution fields for the coarser level
    flecsi::execute<task::rad::const_init<D>, flecsi::default_accelerator>(
      mc, s.Esf(mc), 0.0);
    flecsi::execute<task::rad::const_init<D>, flecsi::default_accelerator>(
      mc, s.Esf(mc, 1), 0.0);

    vcycle<D>(s, index + 1);

    flecsi::execute<task::rad::nlinear_interpolation<D>,
      flecsi::default_accelerator>(mc, mf, s.Esf(mc), s.Errf(mf));
    // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Errf(mf));

    flecsi::execute<task::rad::correction<D>, flecsi::default_accelerator>(
      mf, s.Esf(mf), s.Errf(mf));
    // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Errf(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < s.mg_post; ++i) {
      s.Esf.flip();
      flecsi::execute<task::rad::damped_jacobi<D>, flecsi::default_accelerator>(
        mf, s.Ew(mf), s.Esf(mf), s.Esf(mf, 1), s.Ef(mf), 0.8);
      // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf));
      // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf, 1));
    } // for
  } // if
} // vcycle

template<std::size_t D>
void
fmg(state<D> & s, std::size_t index = 0) {
  auto & mf = *s.mh[index];

  // The scheme requires:
  // 1) Go to a coarser grid, adapt all and repeat this step
  // 2) If in the deeper level, direct solve or do a V-Cycle for a number of

  // iterations 3) Come back up, interpolate, and do a V-Cycle

  // Find current level
  std::size_t level{s.highest_level - index};

  // Deepest level
  if(level == s.lowest_level) {
#ifndef HARD_BENCHMARK_MODE

    flog(warn) << "Deepest level(index):" << level << "(" << index << ")"
               << std::endl;

#endif
    // If in the deepest level, the V-Cycle is already doing a direct solve
    vcycle<D>(s, index);
  }
  else {
#ifndef HARD_BENCHMARK_MODE

    flog(warn) << "cycle level(index): " << level << "(" << index << ")"
               << std::endl;
#endif
    auto & mc = *s.mh[index + 1];

    // Set the RHS and solution field
    flecsi::execute<task::rad::full_weighting<D>, flecsi::default_accelerator>(
      mf, mc, s.Ef(mf), s.Ef(mc));
    flecsi::execute<task::rad::full_weighting<D>, flecsi::default_accelerator>(
      mf, mc, s.Esf(mf), s.Esf(mc));

    // flecsi::execute<task::rad::apply_BC<D>>(mc, s.Ef(mc));
    // flecsi::execute<task::rad::apply_BC<D>>(mc, s.Esf(mc));

    // Set the diffusion coefficient and the stencil (TODO)
    flecsi::execute<task::rad::full_weighting<D>, flecsi::default_accelerator>(
      mf, mc, s.Df_x(mf), s.Df_x(mc));
    // flecsi::execute<task::rad::apply_BC<D>>(mc, s.Df_x(mc));

    flecsi::execute<task::rad::full_weighting<D>, flecsi::default_accelerator>(
      mf, mc, s.Df_y(mf), s.Df_y(mc));
    // flecsi::execute<task::rad::apply_BC<D>>(mc, s.Df_y(mc));

    flecsi::execute<task::rad::full_weighting<D>, flecsi::default_accelerator>(
      mf, mc, s.Df_z(mf), s.Df_z(mc));
    // flecsi::execute<task::rad::apply_BC<D>>(mc, s.Df_z(mc));

    flecsi::execute<task::rad::stencil_init<D>, flecsi::default_accelerator>(
      mc, s.Df_x(mc), s.Df_y(mc), s.Df_z(mc), s.Ew(mc), s.dt(s.gt));

    // Now call solve for one level deeper
    fmg<D>(s, index + 1);

    // Interpolate solution back up (RHS does not change)
    flecsi::execute<task::rad::nlinear_interpolation<D>,
      flecsi::default_accelerator>(mc, mf, s.Esf(mc), s.Esf(mf));
    // flecsi::execute<task::rad::apply_BC<D>>(mf, s.Esf(mf));

    // Do a V-Cycle
    for(std::size_t i{0}; i < s.mg_cycles; ++i) {
      vcycle<D>(s, index);
    } // for
  } // if
} // fmg

} // namespace hard::rad

#endif // HARD_RAD_HH
