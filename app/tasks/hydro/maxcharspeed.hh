
#pragma once

#include "../utils.hh"
#include <cstddef>

namespace hard::tasks::hydro {
template<std::size_t D>
double
update_max_characteristic_speed(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  typename field<vec<D>>::template accessor<ro, na> u_a,
  field<double>::accessor<ro, na> c_a) noexcept {
  using hard::tasks::util::get_mdiota_policy;
  namespace fold = flecsi::exec::fold;

  auto density = m.template mdcolex<is::cells>(r_a);
  auto velocity = m.template mdcolex<is::cells>(u_a);
  auto soundspeed = m.template mdcolex<is::cells>(c_a);

  if constexpr(D == 1) {
    return s.executor().reduceall(
      i, up, (m.template cells<ax::x, dm::quantities>()), fold::max, double) {
      up(velocity(i).norm() + soundspeed(i));
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return s.executor().reduceall(ji, up, mdpolicy_qq, fold::max, double) {
      auto [j, i] = ji;
      up(velocity(i, j).norm() + soundspeed(i, j));
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return s.executor().reduceall(kji, up, mdpolicy_qqq, fold::max, double) {
      auto [k, j, i] = kji;
      up(velocity(i, j, k).norm() + soundspeed(i, j, k));
    };
  }

} //

} // namespace hard::tasks::hydro
