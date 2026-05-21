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
// q[1] = velocity
// q[2] = pressure
// q[3] = specific internal energy
// q[4] = size
// TODO : do we need more field?
template<std::size_t D>
struct extp_q {
  static constexpr std::size_t N_q = D + 4;

  std::array<double, N_q> vals{};

  double & rho() {
    return vals[0];
  }

  spec::vec_view<D> u() {
    return spec::vec_view<D>{vals.data() + 1};
  }

  double & p() {
    return vals[1 + D];
  }
  double & e() {
    return vals[2 + D];
  }
  double & n_comps() {
    return vals[3 + D];
  }

  double rho() const {
    return vals[0];
  }

  vec<D> u() const {
    vec<D> out{};
    for(std::size_t d = 0; d < D; ++d)
      out[d] = vals[1 + d];
    return out;
  }

  double p() const {
    return vals[1 + D];
  }
  double e() const {
    return vals[2 + D];
  }
  double n_comps() const {
    return vals[3 + D];
  }

  double & operator[](std::size_t i) {
    return vals[i];
  }
  const double & operator[](std::size_t i) const {
    return vals[i];
  }

  extp_q & operator=(double x) {
    vals.fill(x);
    return *this;
  }

  std::size_t size() const {
    return N_q;
  }
};

} // namespace hard::tasks::util

#endif // HARD_MODULE_GFM_TASKS_UTIL_HH
