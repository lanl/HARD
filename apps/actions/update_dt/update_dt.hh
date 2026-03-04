#ifndef HARD_APPS_ACTIONS_UPDATE_DT_HH
#define HARD_APPS_ACTIONS_UPDATE_DT_HH

#include "hydro/tasks/maxcharspeed.hh"
#include "hydro/tasks/rhs.hh"

namespace hard::actions {

template<std::size_t D>
void
update_dt(state<D> & s, flecsi::scheduler & sc) {

  auto lmax_f = sc.execute<tasks::hydro::update_max_characteristic_speed<D>>(
    flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.sound_speed(*s.m));

  s.dtmin_ = sc.reduce<tasks::hydro::update_dtmin<D>, flecsi::exec::fold::min>(
    flecsi::exec::on, *s.m, lmax_f);
}

} // namespace hard::actions

#endif // HARD_APPS_ACTIONS_UPDATE_DT_HH
