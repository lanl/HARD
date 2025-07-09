
#pragma once

#include "../../types.hh"
#include <cstddef>

namespace hard::tasks::hydro {

template<typename Limiter>
struct stencil {
  template<typename A, typename T>
  FLECSI_INLINE_TARGET auto operator()(const T & i, const A & acc) {
    return Limiter::reconstruct(
      acc(i - 2), acc(i - 1), acc(i), acc(i + 1), acc(i + 2));
  }

  template<typename A, typename T>
  FLECSI_INLINE_TARGET auto
  operator()(const int & x, const T & i, const T & j, const A & acc) {
    if(x == 0)
      return Limiter::reconstruct(
        acc(i - 2, j), acc(i - 1, j), acc(i, j), acc(i + 1, j), acc(i + 2, j));
    return Limiter::reconstruct(
      acc(i, j - 2), acc(i, j - 1), acc(i, j), acc(i, j + 1), acc(i, j + 2));
  }

  template<typename A, typename T>
  FLECSI_INLINE_TARGET auto operator()(const int & x,
    const T & i,
    const T & j,
    const T & k,
    const A & acc) {
    if(x == 0)
      return Limiter::reconstruct(acc(i - 2, j, k),
        acc(i - 1, j, k),
        acc(i, j, k),
        acc(i + 1, j, k),
        acc(i + 2, j, k));
    if(x == 1)
      return Limiter::reconstruct(acc(i, j - 2, k),
        acc(i, j - 1, k),
        acc(i, j, k),
        acc(i, j + 1, k),
        acc(i, j + 2, k));
    return Limiter::reconstruct(acc(i, j, k - 2),
      acc(i, j, k - 1),
      acc(i, j, k),
      acc(i, j, k + 1),
      acc(i, j, k + 2));
  }
};

namespace utils {

template<typename T>
FLECSI_INLINE_TARGET void
tie(const std::tuple<T, T> & tup, T & a, T & b) noexcept {
  a = std::get<0>(tup);
  b = std::get<1>(tup);
}

} // namespace utils

//
// Perform reconstruction of primitive variables on cell interfaces, calculate
// corresponding conservative variables on faces, and store them into `*Head`
// and `*Tail` variables.
//
//
//    |      U_i     |
//    |       *      |
//     ^HEAD(i)     ^Tail(i)
//    |              |
//   ^Tail(i-1)       ^Head(i+1)
//
//
template<std::size_t Dim, typename Limiter>
void
reconstruct(flecsi::exec::accelerator s,
  std::size_t reconstruction_axis,
  typename mesh<Dim>::template accessor<ro> m,
  // cell-centered primitive varibles
  field<double>::accessor<ro, ro> mass_density_a,
  typename field<vec<Dim>>::template accessor<ro, ro> velocity_a,
  field<double>::accessor<ro, ro> pressure_a,
  field<double>::accessor<ro, ro> specific_internal_energy_a,
  field<double>::accessor<ro, ro> soundspeed_a,
  field<double>::accessor<ro, ro>
#ifdef ENABLE_RADIATION
    radiation_energy_density_a
#endif
  ,
  // reconstructed primitives on faces
  field<double>::accessor<wo, na> rTail_a,
  field<double>::accessor<wo, na> rHead_a,
  typename field<vec<Dim>>::template accessor<wo, na> uTail_a,
  typename field<vec<Dim>>::template accessor<wo, na> uHead_a,
  field<double>::accessor<wo, na> pTail_a,
  field<double>::accessor<wo, na> pHead_a,
  field<double>::accessor<wo, na> eTail_a,
  field<double>::accessor<wo, na> eHead_a,
  field<double>::accessor<wo, na> cTail_a,
  field<double>::accessor<wo, na> cHead_a,
  field<double>::accessor<wo, na>
#ifdef ENABLE_RADIATION
    EradTail_a
#endif
  ,
  field<double>::accessor<wo, na>
#ifdef ENABLE_RADIATION
    EradHead_a
#endif
  ,
  // reconstructed conservatives on faces
  typename field<vec<Dim>>::template accessor<wo, na> ruTail_a,
  typename field<vec<Dim>>::template accessor<wo, na> ruHead_a,
  field<double>::accessor<wo, na> rETail_a,
  field<double>::accessor<wo, na> rEHead_a) noexcept {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto velocity = m.template mdcolex<is::cells>(velocity_a);
  auto pressure = m.template mdcolex<is::cells>(pressure_a);
  auto soundspeed = m.template mdcolex<is::cells>(soundspeed_a);
  auto specific_internal_energy =
    m.template mdcolex<is::cells>(specific_internal_energy_a);

#ifdef ENABLE_RADIATION
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);
#endif
  auto rTail = m.template mdcolex<is::cells>(rTail_a);
  auto rHead = m.template mdcolex<is::cells>(rHead_a);
  auto uTail = m.template mdcolex<is::cells>(uTail_a);
  auto uHead = m.template mdcolex<is::cells>(uHead_a);
  auto pTail = m.template mdcolex<is::cells>(pTail_a);
  auto pHead = m.template mdcolex<is::cells>(pHead_a);
  auto eTail = m.template mdcolex<is::cells>(eTail_a);
  auto eHead = m.template mdcolex<is::cells>(eHead_a);
  auto cTail = m.template mdcolex<is::cells>(cTail_a);
  auto cHead = m.template mdcolex<is::cells>(cHead_a);
#ifdef ENABLE_RADIATION
  auto EradTail = m.template mdcolex<is::cells>(EradTail_a);
  auto EradHead = m.template mdcolex<is::cells>(EradHead_a);
#endif
  auto ruTail = m.template mdcolex<is::cells>(ruTail_a);
  auto ruHead = m.template mdcolex<is::cells>(ruHead_a);
  auto rETail = m.template mdcolex<is::cells>(rETail_a);
  auto rEHead = m.template mdcolex<is::cells>(rEHead_a);

  using hard::tasks::util::get_mdiota_policy;
  using spec::utils::sqr;

  auto ra = reconstruction_axis;

  if constexpr(Dim == 1) {

    s.executor().forall(i, (m.template cells<ax::x, dm::predictor>())) {
      utils::tie(stencil<Limiter>()(i, mass_density), rHead(i), rTail(i));
      utils::tie(stencil<Limiter>()(i, pressure), pHead(i), pTail(i));
      utils::tie(
        stencil<Limiter>()(i, specific_internal_energy), eHead(i), eTail(i));
      utils::tie(stencil<Limiter>()(i, soundspeed), cHead(i), cTail(i));
      utils::tie(stencil<Limiter>()(i, velocity), uHead(i), uTail(i));

      // Compute conservative variables
      ruHead(i) = rHead(i) * uHead(i);
      ruTail(i) = rTail(i) * uTail(i);
      rEHead(i) =
        rHead(i) * eHead(i) + 0.5 * rHead(i) * uHead(i).norm_squared();
      rETail(i) =
        rTail(i) * eTail(i) + 0.5 * rTail(i) * uTail(i).norm_squared();

#ifdef ENABLE_RADIATION
      utils::tie(stencil<Limiter>()(i, radiation_energy_density),
        EradHead(i),
        EradTail(i));
#endif

    }; // forall
  }
  else if constexpr(Dim == 2) {

    auto mdpolicy_pp = get_mdiota_policy(velocity,
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(ji, mdpolicy_pp) {
      auto [j, i] = ji;

      utils::tie(
        stencil<Limiter>()(ra, i, j, mass_density), rHead(i, j), rTail(i, j));
      utils::tie(
        stencil<Limiter>()(ra, i, j, pressure), pHead(i, j), pTail(i, j));
      utils::tie(stencil<Limiter>()(ra, i, j, specific_internal_energy),
        eHead(i, j),
        eTail(i, j));
      utils::tie(
        stencil<Limiter>()(ra, i, j, soundspeed), cHead(i, j), cTail(i, j));
      utils::tie(
        stencil<Limiter>()(ra, i, j, velocity), uHead(i, j), uTail(i, j));

#ifdef ENABLE_RADIATION
      utils::tie(stencil<Limiter>()(ra, i, j, radiation_energy_density),
        EradHead(i, j),
        EradTail(i, j));
#endif

      // Compute conservative variables
      ruHead(i, j) = rHead(i, j) * uHead(i, j);
      ruTail(i, j) = rTail(i, j) * uTail(i, j);
      rEHead(i, j) = rHead(i, j) * eHead(i, j) +
                     0.5 * rHead(i, j) * uHead(i, j).norm_squared();
      rETail(i, j) = rTail(i, j) * eTail(i, j) +
                     0.5 * rTail(i, j) * uTail(i, j).norm_squared();
    }; // forall
  }
  else { // Dim == 3

    auto mdpolicy_ppp = get_mdiota_policy(velocity,
      m.template cells<ax::z, dm::predictor>(),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(kji, mdpolicy_ppp) {
      auto [k, j, i] = kji;

      utils::tie(stencil<Limiter>()(ra, i, j, k, mass_density),
        rHead(i, j, k),
        rTail(i, j, k));
      utils::tie(stencil<Limiter>()(ra, i, j, k, pressure),
        pHead(i, j, k),
        pTail(i, j, k));
      utils::tie(stencil<Limiter>()(ra, i, j, k, specific_internal_energy),
        eHead(i, j, k),
        eTail(i, j, k));
      utils::tie(stencil<Limiter>()(ra, i, j, k, soundspeed),
        cHead(i, j, k),
        cTail(i, j, k));
      utils::tie(stencil<Limiter>()(ra, i, j, k, velocity),
        uHead(i, j, k),
        uTail(i, j, k));

#ifdef ENABLE_RADIATION
      utils::tie(stencil<Limiter>()(ra, i, j, k, radiation_energy_density),
        EradHead(i, j, k),
        EradTail(i, j, k));
#endif

      // Compute conservative variables on faces
      ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);
      ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
      rEHead(i, j, k) = rHead(i, j, k) * eHead(i, j, k) +
                        0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();
      rETail(i, j, k) = eTail(i, j, k) * eTail(i, j, k) +
                        0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
    }; // forall
  }
}

} // namespace hard::tasks::hydro
