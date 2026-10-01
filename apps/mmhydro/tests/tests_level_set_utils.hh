#ifndef TEST_LS_UTILS_HH
#define TEST_LS_UTILS_HH

#include "common/tasks/utils.hh"
#include "spec/labels.hh"
#include "spec/mesh.hh"
#include "spec/types.hh"

using namespace spec;
using namespace common::tasks::utils;

namespace tests::level_set_utils {

template<size_t D>
struct state {
  global::ptr gt;
  typename mesh<D>::ptr m;

  // Velocity field (u)
  static inline const field<vec<D>>::template definition<mesh<D>, is::cells>
    velocity;

  // Time step, time, and normalization period T
  static inline const single<double>::definition<flecsi::topo::global> dt, t,
    period;

  // Level set fields (stage buffers)
  static inline const field<double>::definition<mesh<D>, is::cells> psi;
  static inline const field<double>::definition<mesh<D>, is::cells> psi1, psi2,
    psi3;
  static inline const field<double>::definition<mesh<D>, is::cells> initial_psi;

  // RHS storage for HJ upwind operator: k = L(psi_stage, frozen_material)
  static inline const field<double>::definition<mesh<D>, is::cells> rhs;

  // Material labels (integer IDs)
  static inline const field<int>::definition<mesh<D>, is::cells> material;
  static inline const field<int>::definition<mesh<D>, is::cells>
    frozen_material;
};

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
time_reversal_vel(double t, double period, const std::array<double, D> & x) {

  const double time_factor = std::cos(M_PI * (t / period));

  if constexpr(D == 1) {
    return {2.0 * std::sin(2.0 * M_PI * x[0]) * time_factor};
  }
  else if constexpr(D == 2) {
    return {-2.0 * std::cos(M_PI * x[1]) * std::sin(M_PI * x[0]) *
              std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]) * time_factor,
      2.0 * std::cos(M_PI * x[0]) * std::sin(M_PI * x[0]) *
        std::sin(M_PI * x[1]) * std::sin(M_PI * x[1]) * time_factor};
  }
  else /* D == 3 */ {

    const double sx = std::sin(M_PI * x[0]);
    const double sy = std::sin(M_PI * x[1]);
    const double sz = std::sin(M_PI * x[2]);

    return {2.0 * sx * sx * std::sin(2 * M_PI * x[1]) *
              std::sin(2 * M_PI * x[2]) * time_factor,
      -2.0 * std::sin(2 * M_PI * x[0]) * sy * sy * std::sin(2 * M_PI * x[2]) *
        time_factor,
      -2.0 * std::sin(2 * M_PI * x[0]) * std::sin(2 * M_PI * x[1]) * sz * sz *
        time_factor};
  }
}; // time_reversal_vel

struct copy_test_vel_task : hard::annotation_task_region {
  inline static const std::string name{"copy_test_vel"};
};

template<std::size_t D, int STAGE>
void
copy_test_vel_stage(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<spec::ro> m,
  typename field<vec<D>>::template accessor<spec::wo, spec::na> velocity_a,
  typename single<double>::template accessor<spec::ro> t,
  typename single<double>::template accessor<spec::ro> dt,
  typename single<double>::template accessor<spec::ro> period) noexcept {

  static_assert(STAGE >= 1 && STAGE <= 3, "STAGE must be 1, 2, or 3");

  auto u = m.template mdcolex<is::cells>(velocity_a);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double stage_time = (STAGE == 1)   ? *t
                                : (STAGE == 2) ? (*t + *dt)
                                               : (*t + 0.5 * *dt);
      u(i) = time_reversal_vel<D>(
        stage_time, *period, {m.template center<ax::x>(i)});
    }; // forall
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      const double stage_time = (STAGE == 1)   ? *t
                                : (STAGE == 2) ? (*t + *dt)
                                               : (*t + 0.5 * *dt);

      u(i, j) = time_reversal_vel<D>(stage_time,
        *period,
        {m.template center<ax::x>(i), m.template center<ax::y>(j)});
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      const double stage_time = (STAGE == 1)   ? *t
                                : (STAGE == 2) ? (*t + *dt)
                                               : (*t + 0.5 * *dt);

      u(i, j, k) = time_reversal_vel<D>(stage_time,
        *period,
        {m.template center<ax::x>(i),
          m.template center<ax::y>(j),
          m.template center<ax::z>(k)});
    }; // forall
  }
}; // copy_test_vel_stage

} // namespace tests::level_set_utils

#endif
