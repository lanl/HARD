#ifndef HARD_APPS_MMHYDRO_RECONSTRUCT_DISPATCH_HH
#define HARD_APPS_MMHYDRO_RECONSTRUCT_DISPATCH_HH

#include "gfm/tasks/reconstruct.hh"
#include "spec/labels.hh"
#include "spec/reconstruct.hh"

#include <utility>

namespace hard::actions {

// mmhydro-specific dispatch: selects the reconstructor at runtime and executes
// the material-aware reconstruct_primitives from gfm/tasks/reconstruct.hh.
// Mirrors the hydro dispatch but binds to the gfm overload (which takes the
// material accessor). Kept separate to avoid the overload collision between the
// hydro and gfm reconstruct_primitives, both in hard::tasks::hydro.
template<std::size_t D, typename... Args>
void
dispatch_reconstruct_mm(spec::rc::reconstructor scheme,
  flecsi::scheduler & sc,
  Args &&... args) {
  using namespace spec::reconstructors;
  switch(scheme) {
    case spec::rc::weno5z:
      annotation_rguard<tasks::hydro::reconstruct_primitives_task>(),
      sc.execute<tasks::hydro::reconstruct_primitives<D, weno5z>>(
        flecsi::exec::on, std::forward<Args>(args)...);
      break;
    case spec::rc::ppm4:
      annotation_rguard<tasks::hydro::reconstruct_primitives_task>(),
      sc.execute<tasks::hydro::reconstruct_primitives<D, ppm4>>(
        flecsi::exec::on, std::forward<Args>(args)...);
      break;
    case spec::rc::minmod:
      annotation_rguard<tasks::hydro::reconstruct_primitives_task>(),
      sc.execute<tasks::hydro::reconstruct_primitives<D, minmod>>(
        flecsi::exec::on, std::forward<Args>(args)...);
      break;
    case spec::rc::genminmod:
      annotation_rguard<tasks::hydro::reconstruct_primitives_task>(),
      sc.execute<tasks::hydro::reconstruct_primitives<D, genminmod>>(
        flecsi::exec::on, std::forward<Args>(args)...);
      break;
    case spec::rc::noop:
      annotation_rguard<tasks::hydro::reconstruct_primitives_task>(),
      sc.execute<tasks::hydro::reconstruct_primitives<D, noop>>(
        flecsi::exec::on, std::forward<Args>(args)...);
      break;
  }
}

} // namespace hard::actions

#endif // HARD_APPS_MMHYDRO_RECONSTRUCT_DISPATCH_HH
