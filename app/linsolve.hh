#ifndef LINSOLVE_HH
#define LINSOLVE_HH

#include "flecsi/execution.hh"
#include "flecsi/flog.hh"

#include "rad.hh"
#include "state.hh"

#ifdef USE_FLECSOLVE
#include "flecsolvers.hh"
#endif

namespace hard::rad {

#if defined(USE_FLECSOLVE) && USE_FLECSOLVE

template<std::size_t D>
auto
make_solver(control_policy<state, D> & cp) {

  auto & s = cp.state();
  using namespace flecsolve;

  // Solver parameters
  double temp{0.0};

  // Only pass cp instead?
  solver_parameters<D> params{
    std::ref(cp.state()), std::ref(cp.scheduler()), temp};
  op::core<Ax_op<D>> so(params);
  auto op_handle = op::ref(so);

  precond_parameters<D> pparams{std::ref(s),
    std::ref(cp.scheduler()),
    0,
    s.jacobi_iterations,
    s.nr_vcycles};
  op::core<v_cycle<D>> po(pparams);
  auto prec_handle = op::ref(po);

  auto f = flecsolve::vec::make(s.Ef(*s.m));

  std::size_t iter{0};

  auto slv = flecsolve::bicgstab::solver(
    s.solver_settings, flecsolve::bicgstab::make_work(f))(
    op_handle, prec_handle, [&](auto &, double rnorm) { return false; });
  return slv;
}
#endif

template<std::size_t D>
void
linsolve(control_policy<state, D> & cp) {

  flecsi::scheduler & sc = cp.scheduler();
  auto & s = cp.state();
  auto & mf = *s.m;

  if(s.full_multigrid) {
    fmg<D>(cp);

#if defined(USE_FLECSOLVE) && USE_FLECSOLVE
    sc.execute<task::rad::residual<D>>(
      flecsi::exec::on, mf, s.Ew(mf), s.Esf(mf), s.Ef(mf), s.Resf(mf));
    auto r = flecsolve::vec::make(s.Resf(mf));
    flog(info) << "final res norm radiation: " << r.l2norm().get() << std::endl;
#endif
  }
  else {
#if defined(USE_FLECSOLVE) && USE_FLECSOLVE
    // flecsolve vectors
    auto f = flecsolve::vec::make(s.Ef(mf));
    auto u = flecsolve::vec::make(s.Uf(mf));

    // Normalize rhs
    double f2norm = f.l2norm().get();

    /* Sets up a new linear solver at each iteration (this should not be
     * done).*/
    auto slv = make_solver<D>(cp);
    auto slv_info = slv(f, u);
    auto iters = slv_info.iters;
    auto res_norm_final = slv_info.res_norm_final;

    flog(info) << "final res norm radiation (flecsolve): " << res_norm_final
               << " iter: " << iters << std::endl;

    sc.execute<task::rad::residual<D>>(
      flecsi::exec::on, mf, s.Ew(mf), s.Uf(mf), s.Ef(mf), s.Resf(mf));
    auto r = flecsolve::vec::make(s.Resf(mf));
    flog(info) << "final res norm radiation: " << r.l2norm().get() << std::endl;

#endif
  }
} // linsolve
} // namespace hard::rad

#endif
