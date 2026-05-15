#ifndef HARD_MODULE_GFM_TASKS_UTIL_HH
#define HARD_MODULE_GFM_TASKS_UTIL_HH

#include "../state.hh"

#include <array>
#include <sstream>

namespace hard::tasks::util {

// Define alias extrapolated state. We follow our primitive convention
// q[0] = rho
// q[1...D] = velocity
// q[D+1] = pressure
// q[D+2] = specific internal energy
// TODO : do we need more field?
template<std::size_t D>
using extp_q = std::array<double, D + 3>;

#if 0
template<dm::domain DM, std::size_t D, typename T>
void
copy_scal_field(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<T>::template accessor<ro, ro> from_a,
  typename field<T>::template accessor<wo, wo> to_a) noexcept {

  auto from = m.template mdcolex<is::cells>(from_a);
  auto to = m.template mdcolex<is::cells>(to_a);

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, DM>()) {
      to(i) = from(i);
    } // for
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, DM>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        to(i, j) = from(i, j);
      }
    } // for
  }
  else if constexpr(D == 3) {
    for(auto k : m.template cells<ax::z, DM>()) {
      for(auto j : m.template cells<ax::y, DM>()) {
        for(auto i : m.template cells<ax::x, DM>()) {
          to(i, j, k) = from(i, j, k);
        }
      }
    } // for
  } // if
} // copy_scal_field
#endif

} // namespace hard::tasks::util

#endif // HARD_MODULE_GFM_TASKS_UTIL_HH
