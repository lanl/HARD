#ifndef FLECSOLVERS_HH
#define FLECSOLVERS_HH

#include "flecsi/execution.hh"
#include "flecsi/flog.hh"

#include "state.hh"
#include "tasks/boundaries/boundary.hh"
#include "tasks/rad.hh"
#include "types.hh"

#if defined(USE_FLECSOLVE) && USE_FLECSOLVE
#include "flecsolve/operators/core.hh"
#include "flecsolve/solvers/cg.hh"
#include "flecsolve/solvers/factory.hh"
#include "flecsolve/solvers/gmres.hh"
#include "flecsolve/vectors/topo_view.hh"
#endif

namespace hard::rad {

template<std::size_t D>
struct solver_parameters {
  std::reference_wrapper<state<D>> s;
  std::reference_wrapper<flecsi::scheduler> sc;
  double temp;
};

template<std::size_t D>
struct Ax_op : flecsolve::op::base<solver_parameters<D>> {
  using base = flecsolve::op::base<solver_parameters<D>>;
  using base::params;

  Ax_op(solver_parameters<D> params) : base(std::move(params)) {}

  template<class Domain, class Ranges>
  void apply(const Domain & x, Ranges & y) const {
    flecsi::scheduler & sc = params.sc.get();
    sc.execute<task::rad::Ax<D>>(flecsi::exec::on,
      y.data.topo(),
      std::move(params.s.get().Ew(y.data.topo())),
      y.data.ref(),
      x.data.ref());
    // flecsi::execute<task::rad::apply_radiation_boundary<D>>(flecsi::exec::on,
    //   y.data.topo(),
    //   y.data.ref(),
    //   params.temp); // HARD CODED VALUE for now
  }
};

template<std::size_t D>
struct precond_parameters {
  std::reference_wrapper<state<D>> s;
  std::reference_wrapper<flecsi::scheduler> sc;
  std::size_t index;
  std::size_t nr_vcycles;
  std::size_t jacobi_iterations;
};

template<std::size_t D>
struct v_cycle : flecsolve::op::base<precond_parameters<D>> {
  using base = flecsolve::op::base<precond_parameters<D>>;
  using base::params;

  v_cycle(precond_parameters<D> params) : base(std::move(params)) {}

  template<class Domain, class Range>
  void apply(const Domain & x, Range & y) const {

    flecsi::scheduler & sc = params.sc.get();

    sc.execute<task::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      x.data.ref(),
      params.s.get().Ef_temp(x.data.topo()));

    // Zero solution vector
    sc.execute<task::rad::const_init<D>>(
      flecsi::exec::on, y.data.topo(), params.s.get().Esf(y.data.topo()), 0.0);
    sc.execute<task::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().Esf(y.data.topo(), 1),
      0.0);
    sc.execute<task::rad::const_init<D>>(
      flecsi::exec::on, y.data.topo(), params.s.get().Resf(y.data.topo()), 0.0);

    _vcycle(std::move(params.s), 0);

    sc.execute<task::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().Esf(y.data.topo()),
      y.data.ref());
  }

  void _vcycle(state<D> & s, std::size_t index) const {
    auto & mf = *s.mh[index];
    flecsi::scheduler & sc = params.sc.get();

    // Find current level
    std::size_t level{s.highest_level - index};

    if(level == s.lowest_level) {

      // Direct solve for a single interior point
      for(std::size_t i{0}; i < params.jacobi_iterations; i++) {
        s.Esf.flip();
        // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
        // parallelized
        sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.Ew(mf),
          s.Esf(mf),
          s.Esf(mf, 1),
          s.Ef_temp(mf),
          0.8);
      } // for

      // using namespace flecsolve;
      // solver_parameters<D> params{std::ref(s), std::ref(sc), 0.0};
      // op::core<Ax_op<D>> so(params);

      // auto f = flecsolve::vec::make((s.Ef(mf)));
      // auto u = flecsolve::vec::make((s.Esf(mf)));

      // std::size_t iter{0};
      // auto slv =
      //   flecsolve::cg::solver(s.solver_settings,
      //   flecsolve::cg::make_work(f))(
      //     op::ref(so), op::I, [&](auto &, double rnorm) { return false; });
      // auto info = slv(f, u);
      // flog(info) << "coarse grid res norm " << info.res_norm_final
      //            << " iter: " << info.iters << std::endl;
    }
    else {

      auto & mc = *s.mh[index + 1];
      // Pre Smoothing

      for(std::size_t i{0}; i < s.mg_pre; ++i) {
        s.Esf.flip();
        sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.Ew(mf),
          s.Esf(mf),
          s.Esf(mf, 1),
          s.Ef_temp(mf),
          0.8);
      } // for

      // Set the diffusion coefficient and the stencil (TODO)
      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_x(mf), s.Df_x(mc));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_y(mf), s.Df_y(mc));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_z(mf), s.Df_z(mc));

      sc.execute<task::rad::stencil_init<D>>(flecsi::exec::on,
        mc,
        s.Df_x(mc),
        s.Df_y(mc),
        s.Df_z(mc),
        s.Ew(mc),
        s.dt(*s.gt));

      // Recursive solve
      sc.execute<task::rad::residual<D>>(
        flecsi::exec::on, mf, s.Ew(mf), s.Esf(mf), s.Ef_temp(mf), s.Resf(mf));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Resf(mf), s.Ef_temp(mc));

      // Initialize the solution fields for the coarser level
      sc.execute<task::rad::const_init<D>>(
        flecsi::exec::on, mc, s.Esf(mc), 0.0);
      sc.execute<task::rad::const_init<D>>(
        flecsi::exec::on, mc, s.Esf(mc, 1), 0.0);

      _vcycle(s, index + 1);

      sc.execute<task::rad::cell_centered_interpolation<D>>(
        flecsi::exec::on, mc, mf, s.Esf(mc), s.Errf(mf));

      sc.execute<task::rad::correction<D>>(
        flecsi::exec::on, mf, s.Esf(mf), s.Errf(mf));

      // Post Smoothing
      for(std::size_t i{0}; i < s.mg_post; ++i) {
        s.Esf.flip();
        sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.Ew(mf),
          s.Esf(mf),
          s.Esf(mf, 1),
          s.Ef_temp(mf),
          0.8);
      } // for
    } // if
  }
};

template<std::size_t D>
struct f_mg : flecsolve::op::base<precond_parameters<D>> {
  using base = flecsolve::op::base<precond_parameters<D>>;
  using base::params;

  f_mg(precond_parameters<D> params) : base(std::move(params)) {}

  template<class Domain, class Range>
  void apply(const Domain & x, Range & y) const {

    flecsi::scheduler & sc = params.sc.get();

    // Rhs = r // NO NEED
    sc.execute<task::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      x.data.ref(),
      params.s.get().Ef_temp(x.data.topo()));

    // Zero solution vector
    sc.execute<task::rad::const_init<D>>(
      flecsi::exec::on, y.data.topo(), params.s.get().Esf(y.data.topo()), 0.0);
    sc.execute<task::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().Esf(y.data.topo(), 1),
      0.0);
    sc.execute<task::rad::const_init<D>>(
      flecsi::exec::on, y.data.topo(), params.s.get().Resf(y.data.topo()), 0.0);

    _fmg(std::move(params.s), 0);

    sc.execute<task::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().Esf(y.data.topo()),
      y.data.ref());
  }

  void _vcycle(state<D> & s, std::size_t index) const {
    auto & mf = *s.mh[index];
    flecsi::scheduler & sc = params.sc.get();

    // Find current level
    std::size_t level{s.highest_level - index};

    if(level == s.lowest_level) {

      // Direct solve for a single interior point
      for(std::size_t i{0}; i < params.jacobi_iterations; i++) {
        s.Esf.flip();
        // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
        // parallelized
        sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.Ew(mf),
          s.Esf(mf),
          s.Esf(mf, 1),
          s.Ef_temp(mf),
          0.8);
      } // for

      // using namespace flecsolve;
      // solver_parameters<D> params{std::ref(s), std::ref(sc), 0.0};
      // op::core<Ax_op<D>> so(params);

      // auto f = flecsolve::vec::make((s.Ef(mf)));
      // auto u = flecsolve::vec::make((s.Esf(mf)));

      // std::size_t iter{0};
      // auto slv =
      //   flecsolve::cg::solver(s.solver_settings,
      //   flecsolve::cg::make_work(f))(
      //     op::ref(so), op::I, [&](auto &, double rnorm) { return false; });
      // auto info = slv(f, u);
      // flog(info) << "coarse grid res norm " << info.res_norm_final
      //            << " iter: " << info.iters << std::endl;
    }
    else {

      auto & mc = *s.mh[index + 1];
      // Pre Smoothing

      for(std::size_t i{0}; i < s.mg_pre; ++i) {
        s.Esf.flip();
        sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.Ew(mf),
          s.Esf(mf),
          s.Esf(mf, 1),
          s.Ef_temp(mf),
          0.8);
      } // for

      // Set the diffusion coefficient and the stencil (TODO)
      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_x(mf), s.Df_x(mc));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_y(mf), s.Df_y(mc));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_z(mf), s.Df_z(mc));

      sc.execute<task::rad::stencil_init<D>>(flecsi::exec::on,
        mc,
        s.Df_x(mc),
        s.Df_y(mc),
        s.Df_z(mc),
        s.Ew(mc),
        s.dt(*s.gt));

      // Recursive solve
      sc.execute<task::rad::residual<D>>(
        flecsi::exec::on, mf, s.Ew(mf), s.Esf(mf), s.Ef_temp(mf), s.Resf(mf));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Resf(mf), s.Ef_temp(mc));

      // Initialize the solution fields for the coarser level
      sc.execute<task::rad::const_init<D>>(
        flecsi::exec::on, mc, s.Esf(mc), 0.0);
      sc.execute<task::rad::const_init<D>>(
        flecsi::exec::on, mc, s.Esf(mc, 1), 0.0);

      _vcycle(s, index + 1);

      sc.execute<task::rad::cell_centered_interpolation<D>>(
        flecsi::exec::on, mc, mf, s.Esf(mc), s.Errf(mf));

      sc.execute<task::rad::correction<D>>(
        flecsi::exec::on, mf, s.Esf(mf), s.Errf(mf));

      // Post Smoothing
      for(std::size_t i{0}; i < s.mg_post; ++i) {
        s.Esf.flip();
        sc.execute<task::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.Ew(mf),
          s.Esf(mf),
          s.Esf(mf, 1),
          s.Ef_temp(mf),
          0.8);
      } // for
    } // if
  }

  void _fmg(state<D> & s, std::size_t index) const {
    auto & mf = *s.mh[index];
    flecsi::scheduler & sc = params.sc.get();

    // Find current level
    std::size_t level{s.highest_level - index};

    // Deepest level
    if(level == s.lowest_level) {

      _vcycle(s, index);
    }
    else {
      auto & mc = *s.mh[index + 1];
      // Set the RHS and solution field
      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Ef(mf), s.Ef(mc));
      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Esf(mf), s.Esf(mc));

      // Set the diffusion coefficient and the stencil (TODO)
      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_x(mf), s.Df_x(mc));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_y(mf), s.Df_y(mc));

      sc.execute<task::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.Df_z(mf), s.Df_z(mc));

      sc.execute<task::rad::stencil_init<D>>(flecsi::exec::on,
        mc,
        s.Df_x(mc),
        s.Df_y(mc),
        s.Df_z(mc),
        s.Ew(mc),
        s.dt(*s.gt));

      // Now call solve for one level deeper
      _fmg(s, index + 1);

      // Interpolate solution back up (RHS does not change)
      sc.execute<task::rad::cell_centered_interpolation<D>>(
        flecsi::exec::on, mc, mf, s.Esf(mc), s.Esf(mf));

      // Do a V-Cycle
      for(std::size_t i{0}; i < s.mg_cycles; ++i) {
        _vcycle(s, index);
      } // for
    } // if
  };
};

} // namespace hard::rad

#endif
