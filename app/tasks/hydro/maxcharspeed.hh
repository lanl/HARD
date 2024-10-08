
#pragma once

#include "../utils.hh"
#include <cstddef>

namespace hard::tasks::hydro {
template<std::size_t D>
double
update_max_characteristic_speed(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, na> r_a,
  typename field<vec<D>>::template accessor<ro, na> u_a,
  field<double>::accessor<ro, na> p_a,
  single<double>::accessor<ro> gamma_a) {
  using hard::tasks::util::get_mdiota_policy;
  namespace fold = flecsi::exec::fold;

  auto density = m.template mdcolex<is::cells>(r_a);
  auto velocity = m.template mdcolex<is::cells>(u_a);
  auto pressure = m.template mdcolex<is::cells>(p_a);

  if constexpr(D == 1) {
    return reduceall(i,
      up,
      (m.template cells<ax::x, dm::quantities>()),
      fold::max,
      double,
      "max_char_1d") {
      auto const gamma = *gamma_a;
      const double sound_speed = std::sqrt(gamma * pressure(i) / density(i));
      up(velocity(i).norm() + sound_speed);
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(density,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return reduceall(ji, up, mdpolicy_qq, fold::max, double, "max_char_2d") {
      auto const gamma = *gamma_a;
      auto [j, i] = ji;
      const double sound_speed =
        std::sqrt(gamma * pressure(i, j) / density(i, j));
      up(velocity(i, j).norm() + sound_speed);
    };
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(density,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    return reduceall(kji, up, mdpolicy_qqq, fold::max, double, "max_char_3d") {
      auto const gamma = *gamma_a;
      auto [k, j, i] = kji;

      const double sound_speed =
        std::sqrt(gamma * pressure(i, j, k) / density(i, j, k));
      up(velocity(i, j, k).norm() + sound_speed);
    };
  }

} //

} // namespace hard::tasks::hydro
