
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
reconstruct_primitives(flecsi::exec::accelerator s,
  std::size_t reconstruction_axis,
  typename mesh<Dim>::template accessor<ro> m,
  // cell-centered primitive varibles
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    typename faces<Dim>::accessor<wo, na>>> cons_faces_a,
  std::vector<std::tuple<typename field<vec<Dim>>::accessor<ro, na>,
    typename faces_vec<Dim>::accessor<wo, na>>> cons_faces_vec_a) noexcept {

  using hard::tasks::util::get_mdiota_policy;
  using spec::utils::sqr;

  auto ra = reconstruction_axis;

  if constexpr(Dim == 1) {

    for(auto cf : cons_faces_a) {
      auto [cons_a, face_a] = cf;
      auto [h, t] = faces<Dim>::mdcolex(m, face_a);
      auto cons = m.template mdcolex<is::cells>(cons_a);
      s.executor().forall(i, (m.template cells<ax::x, dm::predictor>())) {
        utils::tie(stencil<Limiter>()(i, cons), h(i), t(i));
      };
    }

    for(auto cf : cons_faces_vec_a) {
      auto [cons_a, face_a] = cf;
      auto [h, t] = faces<Dim>::mdcolex(m, face_a);
      auto cons = m.template mdcolex<is::cells>(cons_a);
      s.executor().forall(i, (m.template cells<ax::x, dm::predictor>())) {
        utils::tie(stencil<Limiter>()(i, cons), h(i), t(i));
      };
    }
  }
  else if constexpr(Dim == 2) {

    auto mdpolicy_pp = get_mdiota_policy(
      m.template mdcolex<is::cells>(std::get<0>(cons_faces_a.front())),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    for(auto cf : cons_faces_a) {
      auto [cons_a, face_a] = cf;
      auto [h, t] = faces<Dim>::mdcolex(m, face_a);
      auto cons = m.template mdcolex<is::cells>(cons_a);
      s.executor().forall(ji, mdpolicy_pp) {
        auto [j, i] = ji;
        utils::tie(stencil<Limiter>()(ra, i, j, cons), h(i, j), t(i, j));
      };
    };
    for(auto cf : cons_faces_vec_a) {
      auto [cons_a, face_a] = cf;
      auto [h, t] = faces<Dim>::mdcolex(m, face_a);
      auto cons = m.template mdcolex<is::cells>(cons_a);
      s.executor().forall(ji, mdpolicy_pp) {
        auto [j, i] = ji;
        utils::tie(stencil<Limiter>()(ra, i, j, cons), h(i, j), t(i, j));
      };
    };
  }
  else { // Dim == 3

    auto mdpolicy_ppp = get_mdiota_policy(
      m.template mdcolex<is::cells>(std::get<0>(cons_faces_a.front())),
      m.template cells<ax::z, dm::predictor>(),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    for(auto cf : cons_faces_a) {
      auto [cons_a, face_a] = cf;
      auto [h, t] = faces<Dim>::mdcolex(m, face_a);
      auto cons = m.template mdcolex<is::cells>(cons_a);
      s.executor().forall(kji, mdpolicy_ppp) {
        auto [k, j, i] = kji;
        utils::tie(
          stencil<Limiter>()(ra, i, j, k, cons), h(i, j, k), t(i, j, k));
      };
    }

    for(auto cf : cons_faces_vec_a) {
      auto [cons_a, face_a] = cf;
      auto [h, t] = faces<Dim>::mdcolex(m, face_a);
      auto cons = m.template mdcolex<is::cells>(cons_a);
      s.executor().forall(kji, mdpolicy_ppp) {
        auto [k, j, i] = kji;
        utils::tie(
          stencil<Limiter>()(ra, i, j, k, cons), h(i, j, k), t(i, j, k));
      };
    }
  }
}

template<std::size_t Dim>
void
reconstruct_conservatives(flecsi::exec::accelerator s,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces<Dim>::accessor<ro, na> rFace_a,
  typename faces_vec<Dim>::accessor<ro, na> uFace_a,
  typename faces<Dim>::accessor<ro, na> eFace_a,
  // reconstructed conservatives on faces
  typename faces_vec<Dim>::accessor<wo, na> ruFace_a,
  typename faces<Dim>::accessor<wo, na> rEFace_a) noexcept {

  auto [rHead, rTail] = faces<Dim>::mdcolex(m, rFace_a);
  auto [uHead, uTail] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [eHead, eTail] = faces<Dim>::mdcolex(m, eFace_a);
  auto [ruHead, ruTail] = faces_vec<Dim>::mdcolex(m, ruFace_a);
  auto [rEHead, rETail] = faces<Dim>::mdcolex(m, rEFace_a);

  using hard::tasks::util::get_mdiota_policy;
  using spec::utils::sqr;

  if constexpr(Dim == 1) {

    s.executor().forall(i, (m.template cells<ax::x, dm::predictor>())) {

      // Compute conservative variables
      ruHead(i) = rHead(i) * uHead(i);
      ruTail(i) = rTail(i) * uTail(i);
      rEHead(i) =
        rHead(i) * eHead(i) + 0.5 * rHead(i) * uHead(i).norm_squared();
      rETail(i) =
        rTail(i) * eTail(i) + 0.5 * rTail(i) * uTail(i).norm_squared();

    }; // forall
  }
  else if constexpr(Dim == 2) {

    auto mdpolicy_pp = get_mdiota_policy(rHead,
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(ji, mdpolicy_pp) {
      auto [j, i] = ji;

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

    auto mdpolicy_ppp = get_mdiota_policy(rHead,
      m.template cells<ax::z, dm::predictor>(),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(kji, mdpolicy_ppp) {
      auto [k, j, i] = kji;

      // Compute conservative variables on faces
      ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);
      ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
      rEHead(i, j, k) = rHead(i, j, k) * eHead(i, j, k) +
                        0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();
      rETail(i, j, k) = rTail(i, j, k) * eTail(i, j, k) +
                        0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
    }; // forall
  }
}

} // namespace hard::tasks::hydro
