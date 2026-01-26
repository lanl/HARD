#ifndef TEST_LS_UTILS_HH
#define TEST_LS_UTILS_HH

#include "../tasks/utils.hh"

using namespace spec;
using namespace hard::tasks::util;

namespace tests::level_set_utils {

template<size_t D>
struct state {
  global::ptr gt;
  typename mesh<D>::ptr m;

  // Velocity field (u)
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    velocity;

  // Time step and time
  static inline const single<double>::definition<flecsi::topo::global> dt, t;

  // Level set fields (stage buffers)
  static inline const field<double>::definition<mesh<D>, is::cells> psi;
  static inline const field<double>::definition<mesh<D>, is::cells> psi1, psi2,
    psi3;

  // RHS storage for HJ upwind operator: k = L(psi_stage, frozen_material)
  static inline const field<double>::definition<mesh<D>, is::cells> rhs;

  // Material labels (integer IDs)
  static inline const field<int>::definition<mesh<D>, is::cells> material;
  static inline const field<int>::definition<mesh<D>, is::cells>
    frozen_material;

  // For unit test comparison
  static inline const field<int>::definition<mesh<D>, is::cells>
    initial_material;
};

template<std::size_t D>
vec<D>
time_reversal_vel(double t, const std::array<double, D> & x) {

  if constexpr(D == 1) {
    return {2.0 * std::sin(2.0 * M_PI * x[0]) * std::cos(M_PI * t / 2.0)};
  }
  else if constexpr(D == 2) {
    return {-2.0 * std::cos(M_PI * x[1]) * std::sin(M_PI * x[0]) *
              std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]) *
              std::cos(M_PI * t / 2.0),
      2.0 * std::cos(M_PI * x[0]) * std::sin(M_PI * x[0]) *
        std::sin(M_PI * x[1]) * std::sin(M_PI * x[1]) *
        std::cos(M_PI * t / 2.0)};
  }
  else /* D == 3 */ {

    const double sx = std::sin(M_PI * x[0]);
    const double sy = std::sin(M_PI * x[1]);
    const double sz = std::sin(M_PI * x[2]);

    return {2.0 * sx * sx * std::sin(2 * M_PI * x[1]) *
              std::sin(2 * M_PI * x[2]) * std::cos(M_PI * t),
      -2.0 * std::sin(2 * M_PI * x[0]) * sy * sy * std::sin(2 * M_PI * x[2]) *
        std::cos(M_PI * t),
      -2.0 * std::sin(2 * M_PI * x[0]) * std::sin(2 * M_PI * x[1]) * sz * sz *
        std::cos(M_PI * t)};
  }
}; // time_reversal_vel

template<std::size_t D>
void
copy_test_vel(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<vec<D>>::template accessor<wo, ro> velocity_a,
  typename single<double>::template accessor<ro> t) noexcept {

  auto u = m.template mdcolex<is::cells>(velocity_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      u(i) = time_reversal_vel<D>(t, {m.template center<ax::x>(i)});
    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;

      u(i, j) = time_reversal_vel<D>(
        t, {m.template center<ax::x>(i), m.template center<ax::y>(j)});
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;

      u(i, j, k) = time_reversal_vel<D>(t,
        {m.template center<ax::x>(i),
          m.template center<ax::y>(j),
          m.template center<ax::z>(k)});
    }; // forall
  }
}; // copy_test_vel

} // namespace tests::level_set_utils

#endif
