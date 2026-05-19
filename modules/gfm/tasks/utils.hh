#ifndef HARD_MODULE_GFM_TASKS_UTIL_HH
#define HARD_MODULE_GFM_TASKS_UTIL_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/limiter.hh"
#include "../../spec/mesh.hh"
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

template<std::size_t D>
struct extp_q_ids {
  static constexpr std::size_t rho = 0;
  static constexpr std::size_t vel = 1;
  static constexpr std::size_t pressure = 1 + D;
  static constexpr std::size_t sie = 2 + D;
  static constexpr std::size_t size = D + 3;
};

} // namespace hard::tasks::util

#endif // HARD_MODULE_GFM_TASKS_UTIL_HH
