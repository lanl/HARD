#ifndef HARD_MODULE_GFM_TASKS_UTIL_HH
#define HARD_MODULE_GFM_TASKS_UTIL_HH

#include "../../common/tasks/utils.hh"
#include "../../hydro/state.hh"
#include "../../spec/labels.hh"
#include "../../spec/limiter.hh"
#include "../../spec/mesh.hh"

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

  std::array<double, N_q> vals;

  FLECSI_INLINE_TARGET double & rho() {
    return vals[0];
  }

  FLECSI_INLINE_TARGET spec::vec_view<D> u() {
    return spec::vec_view<D>{vals.data() + 1};
  }

  FLECSI_INLINE_TARGET double & p() {
    return vals[1 + D];
  }
  FLECSI_INLINE_TARGET double & e() {
    return vals[2 + D];
  }
  FLECSI_INLINE_TARGET double & n_comps() {
    return vals[3 + D];
  }

  FLECSI_INLINE_TARGET double rho() const {
    return vals[0];
  }

  FLECSI_INLINE_TARGET vec<D> u() const {
    vec<D> out{};
    for(std::size_t d = 0; d < D; ++d)
      out[d] = vals[1 + d];
    return out;
  }

  FLECSI_INLINE_TARGET double p() const {
    return vals[1 + D];
  }
  FLECSI_INLINE_TARGET double e() const {
    return vals[2 + D];
  }
  FLECSI_INLINE_TARGET double n_comps() const {
    return vals[3 + D];
  }

  FLECSI_INLINE_TARGET double & operator[](std::size_t i) {
    return vals[i];
  }
  FLECSI_INLINE_TARGET const double & operator[](std::size_t i) const {
    return vals[i];
  }

  FLECSI_INLINE_TARGET extp_q & operator=(double x) {
    vals.fill(x);
    return *this;
  }

  FLECSI_INLINE_TARGET std::size_t size() const {
    return N_q;
  }
};

static_assert(std::is_trivially_default_constructible_v<extp_q<1>>);
static_assert(std::is_trivially_copyable_v<extp_q<1>>);

// ------
template<std::size_t D>
void
material_to_scalar_output(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<int>::template accessor<ro, ro> material_a,
  typename field<double>::template accessor<wo, na> material_out_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto material = m.template mdcolex<is::cells>(material_a);
  auto material_out = m.template mdcolex<is::cells>(material_out_a);

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      material_out(i) = static_cast<double>(material(i));
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(material_out,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      material_out(i, j) = static_cast<double>(material(i, j));
    };
  }
  else {
    auto pol = get_mdiota_policy(material_out,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      material_out(i, j, k) = static_cast<double>(material(i, j, k));
    };
  }
}
// -------

} // namespace hard::tasks::util

#endif // HARD_MODULE_GFM_TASKS_UTIL_HH
