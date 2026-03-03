#ifndef HARD_MODULES_RAD_FLECSOLVERS_HH
#define HARD_MODULES_RAD_FLECSOLVERS_HH

#include "flecsi/execution.hh"
#include "flecsi/flog.hh"

#include "../modules/rad/tasks/rad.hh"

#include "../modules/spec/tasks/boundaries/boundary.hh"

#include "flecsolve/operators/core.hh"
#include "flecsolve/solvers/cg.hh"
#include "flecsolve/solvers/factory.hh"
#include "flecsolve/solvers/gmres.hh"
#include "flecsolve/vectors/topo_view.hh"

namespace hard {

template<std::size_t D>
struct solver_parameters {
  std::reference_wrapper<state<D>> s;
  std::reference_wrapper<flecsi::scheduler> sc;
  double temp;
};

template<std::size_t D>
struct operator_t : flecsolve::op::base<solver_parameters<D>> {
  using base = flecsolve::op::base<solver_parameters<D>>;
  using base::params;

  operator_t(solver_parameters<D> params) : base(std::move(params)) {}

  template<class Domain, class Ranges>
  void apply(const Domain & x, Ranges & y) const {
    flecsi::scheduler & sc = params.sc.get();
    sc.execute<tasks::rad::apply_operator<D>>(flecsi::exec::on,
      y.data.topo(),
      std::move(params.s.get().rad.mgr.Ew(y.data.topo())),
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

    sc.execute<tasks::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      x.data.ref(),
      params.s.get().rad.mgr.Ef_temp(x.data.topo()));

    // Zero solution vector
    sc.execute<tasks::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().rad.mgr.Esf(y.data.topo()),
      0.0);
    sc.execute<tasks::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().rad.mgr.Esf(y.data.topo(), 1),
      0.0);
    sc.execute<tasks::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().rad.mgr.Resf(y.data.topo()),
      0.0);

    _vcycle(std::move(params.s), 0);

    sc.execute<tasks::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().rad.mgr.Esf(y.data.topo()),
      y.data.ref());
  }

  void _vcycle(state<D> & s, std::size_t index) const {
    auto & mf = *s.mh[index];
    flecsi::scheduler & sc = params.sc.get();

    // Find current level
    std::size_t level{s.min_highest_level - index};

    if(level == s.lowest_level) {

      for(std::size_t i{0}; i < params.jacobi_iterations; i++) {
        s.rad.mgr.Esf.flip();
        // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
        // parallelized
        sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.rad.mgr.Ew(mf),
          s.rad.mgr.Esf(mf),
          s.rad.mgr.Esf(mf, 1),
          s.rad.mgr.Ef_temp(mf),
          0.8);
      } // for

      // using namespace flecsolve;
      // solver_parameters<D> params{std::ref(s), std::ref(sc), 0.0};
      // op::core<operator_t<D>> so(params);

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

      for(std::size_t i{0}; i < s.rad.mgr.mg_pre; ++i) {
        s.rad.mgr.Esf.flip();
        sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.rad.mgr.Ew(mf),
          s.rad.mgr.Esf(mf),
          s.rad.mgr.Esf(mf, 1),
          s.rad.mgr.Ef_temp(mf),
          0.8);
      } // for

      // Set the diffusion coefficient and the stencil (TODO)
      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.rad.mgr.Df_x(mf), s.rad.mgr.Df_x(mc));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.rad.mgr.Df_y(mf), s.rad.mgr.Df_y(mc));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.rad.mgr.Df_z(mf), s.rad.mgr.Df_z(mc));

      sc.execute<tasks::rad::stencil_init<D>>(flecsi::exec::on,
        mc,
        s.rad.mgr.Df_x(mc),
        s.rad.mgr.Df_y(mc),
        s.rad.mgr.Df_z(mc),
        s.rad.mgr.Ew(mc),
        s.dt(*s.gt));

      // Recursive solve
      sc.execute<tasks::rad::residual<D>>(flecsi::exec::on,
        mf,
        s.rad.mgr.Ew(mf),
        s.rad.mgr.Esf(mf),
        s.rad.mgr.Ef_temp(mf),
        s.rad.mgr.Resf(mf));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.rad.mgr.Resf(mf), s.rad.mgr.Ef_temp(mc));

      // Initialize the solution fields for the coarser level
      sc.execute<tasks::rad::const_init<D>>(
        flecsi::exec::on, mc, s.rad.mgr.Esf(mc), 0.0);
      sc.execute<tasks::rad::const_init<D>>(
        flecsi::exec::on, mc, s.rad.mgr.Esf(mc, 1), 0.0);

      _vcycle(s, index + 1);

      sc.execute<tasks::rad::cell_centered_interpolation<D>>(
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
          s.rad.mgr.Ef_temp(mf),
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
    sc.execute<tasks::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      x.data.ref(),
      params.s.get().Ef_temp(x.data.topo()));

    // Zero solution vector
    sc.execute<tasks::rad::const_init<D>>(
      flecsi::exec::on, y.data.topo(), params.s.get().Esf(y.data.topo()), 0.0);
    sc.execute<tasks::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().mgr.Esf(y.data.topo(), 1),
      0.0);
    sc.execute<tasks::rad::const_init<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().mgr.Resf(y.data.topo()),
      0.0);

    _fmg(std::move(params.s), 0);

    sc.execute<tasks::rad::copy_field<D>>(flecsi::exec::on,
      y.data.topo(),
      params.s.get().mgr.Esf(y.data.topo()),
      y.data.ref());
  }

  void _vcycle(state<D> & s, std::size_t index) const {
    auto & mf = *s.mh[index];
    flecsi::scheduler & sc = params.sc.get();

    // Find current level
    std::size_t level{s.mgr.min_highest_level - index};

    if(level == s.mgr.lowest_level) {

      for(std::size_t i{0}; i < params.mgr.jacobi_iterations; i++) {
        s.Esf.flip();
        // NOTE: We are defaulting to damped_jacobi until gauss-seidel is
        // parallelized
        sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.mgr.Ew(mf),
          s.mgr.Esf(mf),
          s.mgr.Esf(mf, 1),
          s.mgr.Ef_temp(mf),
          0.8);
      } // for

      // using namespace flecsolve;
      // solver_parameters<D> params{std::ref(s), std::ref(sc), 0.0};
      // op::core<operator_t<D>> so(params);

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

      for(std::size_t i{0}; i < s.mgr.mg_pre; ++i) {
        s.mgr.Esf.flip();
        sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.mgr.Ew(mf),
          s.mgr.Esf(mf),
          s.mgr.Esf(mf, 1),
          s.mgr.Ef_temp(mf),
          0.8);
      } // for

      // Set the diffusion coefficient and the stencil (TODO)
      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Df_x(mf), s.mgr.Df_x(mc));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Df_y(mf), s.mgr.Df_y(mc));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Df_z(mf), s.mgr.Df_z(mc));

      sc.execute<tasks::rad::stencil_init<D>>(flecsi::exec::on,
        mc,
        s.mgr.Df_x(mc),
        s.mgr.Df_y(mc),
        s.mgr.Df_z(mc),
        s.mgr.Ew(mc),
        s.dt(*s.gt));

      // Recursive solve
      sc.execute<tasks::rad::residual<D>>(flecsi::exec::on,
        mf,
        s.mgr.Ew(mf),
        s.mgr.Esf(mf),
        s.mgr.Ef_temp(mf),
        s.mgr.Resf(mf));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Resf(mf), s.mgr.Ef_temp(mc));

      // Initialize the solution fields for the coarser level
      sc.execute<tasks::rad::const_init<D>>(
        flecsi::exec::on, mc, s.mgr.Esf(mc), 0.0);
      sc.execute<tasks::rad::const_init<D>>(
        flecsi::exec::on, mc, s.mgr.Esf(mc, 1), 0.0);

      _vcycle(s, index + 1);

      sc.execute<tasks::rad::cell_centered_interpolation<D>>(
        flecsi::exec::on, mc, mf, s.mgr.Esf(mc), s.mgr.Errf(mf));

      sc.execute<tasks::rad::correction<D>>(
        flecsi::exec::on, mf, s.mgr.Esf(mf), s.Errf(mf));

      // Post Smoothing
      for(std::size_t i{0}; i < s.mg_post; ++i) {
        s.mgr.Esf.flip();
        sc.execute<tasks::rad::damped_jacobi<D>>(flecsi::exec::on,
          mf,
          s.mgr.Ew(mf),
          s.mgr.Esf(mf),
          s.mgr.Esf(mf, 1),
          s.mgr.Ef_temp(mf),
          0.8);
      } // for
    } // if
  }

  void _fmg(state<D> & s, std::size_t index) const {
    auto & mf = *s.mh[index];
    flecsi::scheduler & sc = params.sc.get();

    // Find current level
    std::size_t level{s.mgr.min_highest_level - index};

    // Deepest level
    if(level == s.mgr.lowest_level) {

      _vcycle(s, index);
    }
    else {
      auto & mc = *s.mh[index + 1];
      // Set the RHS and solution field
      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Ef(mf), s.mgr.Ef(mc));
      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Esf(mf), s.mgr.Esf(mc));

      // Set the diffusion coefficient and the stencil (TODO)
      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Df_x(mf), s.mgr.Df_x(mc));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Df_y(mf), s.mgr.Df_y(mc));

      sc.execute<tasks::rad::cell_centered_weighting<D>>(
        flecsi::exec::on, mf, mc, s.mgr.Df_z(mf), s.mgr.Df_z(mc));

      sc.execute<tasks::rad::stencil_init<D>>(flecsi::exec::on,
        mc,
        s.mgr.Df_x(mc),
        s.mgr.Df_y(mc),
        s.mgr.Df_z(mc),
        s.mgr.Ew(mc),
        s.dt(*s.gt));

      // Now call solve for one level deeper
      _fmg(s, index + 1);

      // Interpolate solution back up (RHS does not change)
      sc.execute<tasks::rad::cell_centered_interpolation<D>>(
        flecsi::exec::on, mc, mf, s.mgr.Esf(mc), s.mgr.Esf(mf));

      // Do a V-Cycle
      for(std::size_t i{0}; i < s.mgr.mg_cycles; ++i) {
        _vcycle(s, index);
      } // for
    } // if
  };
};

} // namespace hard

#endif
