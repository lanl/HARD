#ifndef LINSOLVE_HH
#define LINSOLVE_HH

#include "flecsi/execution.hh"
#include "flecsi/flog.hh"

#include "flecsolvers.hh"

namespace hard::rad {

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
  op::core<operator_t<D>> so(params);
  auto op_handle = op::ref(so);

  precond_parameters<D> pparams{std::ref(s),
    std::ref(cp.scheduler()),
    0,
    s.mgr.jacobi_iterations,
    s.mgr.nr_vcycles};
  op::core<v_cycle<D>> po(pparams);
  auto prec_handle = op::ref(po);

  auto f = flecsolve::vec::make(s.mgr.Ef(*s.m));

  auto slv = flecsolve::bicgstab::solver(
    s.mgr.solver_settings, flecsolve::bicgstab::make_work(f))(
    op_handle, prec_handle, [&](auto &, double) { return false; });
  return slv;
}

template<std::size_t D>
void
linsolve(control_policy<state, D> & cp) {

  flecsi::scheduler & sc = cp.scheduler();
  auto & s = cp.state();
  auto & mf = *s.m;

  if(s.mgr.full_multigrid) {
    fmg<D>(cp);

    sc.execute<tasks::rad::residual<D>>(flecsi::exec::on,
      mf,
      s.mgr.Ew(mf),
      s.mgr.Esf(mf),
      s.mgr.Ef(mf),
      s.mgr.Resf(mf));
    auto r = flecsolve::vec::make(s.mgr.Resf(mf));
    flog(info) << "final res norm radiation: " << r.l2norm().get() << std::endl;
  }
  else {
    // flecsolve vectors
    auto f = flecsolve::vec::make(s.mgr.Ef(mf));
    auto u = flecsolve::vec::make(s.mgr.Uf(mf));

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
      s.mgr.Ew(mf),
      s.mgr.Uf(mf),
      s.mgr.Ef(mf),
      s.mgr.Resf(mf));
    auto r = flecsolve::vec::make(s.mgr.Resf(mf));
    flog(info) << "final res norm radiation: " << r.l2norm().get() << std::endl;
  }
} // linsolve
} // namespace hard::rad

#endif
