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
  typename field<int>::template accessor<ro, na> material_a,
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

template<std::size_t D>
void
extp_rho_to_scalar_output(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro>,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<ro, na> q_a,
  typename field<double>::template accessor<wo, na> rho_out_a) noexcept {

  s.executor().forall(i, flecsi::util::iota_view({}, q_a.span().size())) {
    rho_out_a(i) = q_a(i).rho();
  };
}

template<std::size_t D>
void
vec_components_to_scalar_output(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<vec<D>>::template accessor<ro, na> v_a,
  typename field<double>::template accessor<wo, na> x_out_a,
  typename field<double>::template accessor<wo, na> y_out_a,
  typename field<double>::template accessor<wo, na> z_out_a,
  typename field<double>::template accessor<wo, na> mag_out_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  auto v = m.template mdcolex<is::cells>(v_a);
  auto x_out = m.template mdcolex<is::cells>(x_out_a);
  auto y_out = m.template mdcolex<is::cells>(y_out_a);
  auto z_out = m.template mdcolex<is::cells>(z_out_a);
  auto mag_out = m.template mdcolex<is::cells>(mag_out_a);

  if constexpr(D == 1) {
    s.executor().forall(i, m.template cells<ax::x, dm::quantities>()) {
      x_out(i) = v(i).x();
      y_out(i) = 0.0;
      z_out(i) = 0.0;
      mag_out(i) = v(i).norm();
    };
  }
  else if constexpr(D == 2) {
    auto pol = get_mdiota_policy(v,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, pol) {
      auto [j, i] = ji;
      x_out(i, j) = v(i, j).x();
      y_out(i, j) = v(i, j).y();
      z_out(i, j) = 0.0;
      mag_out(i, j) = v(i, j).norm();
    };
  }
  else {
    auto pol = get_mdiota_policy(v,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, pol) {
      auto [k, j, i] = kji;
      x_out(i, j, k) = v(i, j, k).x();
      y_out(i, j, k) = v(i, j, k).y();
      z_out(i, j, k) = v(i, j, k).z();
      mag_out(i, j, k) = v(i, j, k).norm();
    };
  }
}
// -------

} // namespace hard::tasks::util

#endif // HARD_MODULE_GFM_TASKS_UTIL_HH
