#ifndef HARD_APPS_ACTIONS_RHS_HH
#define HARD_APPS_ACTIONS_RHS_HH

#include "hydro/tasks/init.hh"
#include "hydro/tasks/rhs.hh"

namespace hard::actions {

template<std::size_t D>
using field_double =
  flecsi::field<double>::Reference<spec::mesh<D>, spec::is::cells>;

template<std::size_t D>
void
rhs(state<D> & s,
  flecsi::scheduler & sc,
  // conservative, dt1, dt2, n
  std::vector<std::
      tuple<field_double<D>, field_double<D>, field_double<D>, field_double<D>>>
    v_f = {}) {

  sc.execute<tasks::init::compute_dt_weighted>(flecsi::exec::on,
    s.dt(*s.gt),
    s.dt_weighted(*s.gt),
    hard::time_stepper::time_stepper_gamma);

  {

    auto scalar_v = std::vector{// dt1
      s.rk_dt1.mass_density()(*s.m),
      s.rk_dt1.total_energy_density()(*s.m),
      // dt2
      s.rk_dt2.mass_density()(*s.m),
      s.rk_dt2.total_energy_density()(*s.m)};

    for(auto & v : v_f) {
      scalar_v.push_back(std::get<1>(v));
      scalar_v.push_back(std::get<2>(v));
    }

    // Set all dU_dt temporaries to zero before adding time derivative terms
    sc.execute<tasks::hydro::set_dudt_to_zero<D>>(flecsi::exec::on,
      scalar_v,
      std::vector{// dt1
        s.rk_dt1.momentum_energy_density()(*s.m),
        // dt2
        s.rk_dt2.momentum_energy_density()(*s.m)});
  }

  {

    auto scalar_v = std::vector{std::make_tuple(s.cons.hydro.mass_density(*s.m),
                                  s.rk_n.mass_density()(*s.m)),
      std::make_tuple(s.cons.hydro.total_energy_density(*s.m),
        s.rk_n.total_energy_density()(*s.m))};

    for(auto & v : v_f) {
      scalar_v.push_back(std::make_tuple(std::get<0>(v), std::get<3>(v)));
    }

    // Store the current state of evolved variables (U^n) before performing a
    // time step
    // clang-format off
  sc.execute<tasks::hydro::store_current_state<D>>(flecsi::exec::on,
    scalar_v,
    std::vector{
      std::make_tuple(
        s.cons.hydro.momentum_density(*s.m),
        s.rk_n.momentum_energy_density()(*s.m))});
    // clang-format on
  }
}

} // namespace hard::actions

#endif // HARD_APPS_ACTIONS_RHS_HH
