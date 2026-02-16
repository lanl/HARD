#ifndef HARD_MODULE_HYDRO_RECONSTRUCT_HH
#define HARD_MODULE_HYDRO_RECONSTRUCT_HH

#include "types.hh"
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
// corresponding conservative variables on faces, and store them into `*right`
// and `*left` variables.
//
//
//    |      U_i     |
//    |       *      |
//     ^right(i)     ^left(i)
//    |              |
//   ^left(i-1)       ^right(i+1)
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

  auto [rRight, rLeft] = faces<Dim>::mdcolex(m, rFace_a);
  auto [uRight, uLeft] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [eRight, eLeft] = faces<Dim>::mdcolex(m, eFace_a);
  auto [ruRight, ruLeft] = faces_vec<Dim>::mdcolex(m, ruFace_a);
  auto [rERight, rELeft] = faces<Dim>::mdcolex(m, rEFace_a);

  using hard::tasks::util::get_mdiota_policy;
  using spec::utils::sqr;

  if constexpr(Dim == 1) {

    s.executor().forall(i, (m.template cells<ax::x, dm::predictor>())) {

      // Compute conservative variables
      ruRight(i) = rRight(i) * uRight(i);
      ruLeft(i) = rLeft(i) * uLeft(i);
      rERight(i) =
        rRight(i) * eRight(i) + 0.5 * rRight(i) * uRight(i).norm_squared();
      rELeft(i) =
        rLeft(i) * eLeft(i) + 0.5 * rLeft(i) * uLeft(i).norm_squared();

    }; // forall
  }
  else if constexpr(Dim == 2) {

    auto mdpolicy_pp = get_mdiota_policy(rRight,
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(ji, mdpolicy_pp) {
      auto [j, i] = ji;

      // Compute conservative variables
      ruRight(i, j) = rRight(i, j) * uRight(i, j);
      ruLeft(i, j) = rLeft(i, j) * uLeft(i, j);
      rERight(i, j) = rRight(i, j) * eRight(i, j) +
                      0.5 * rRight(i, j) * uRight(i, j).norm_squared();
      rELeft(i, j) = rLeft(i, j) * eLeft(i, j) +
                     0.5 * rLeft(i, j) * uLeft(i, j).norm_squared();
    }; // forall
  }
  else { // Dim == 3

    auto mdpolicy_ppp = get_mdiota_policy(rRight,
      m.template cells<ax::z, dm::predictor>(),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(kji, mdpolicy_ppp) {
      auto [k, j, i] = kji;

      // Compute conservative variables on faces
      ruRight(i, j, k) = rRight(i, j, k) * uRight(i, j, k);
      ruLeft(i, j, k) = rLeft(i, j, k) * uLeft(i, j, k);
      rERight(i, j, k) = rRight(i, j, k) * eRight(i, j, k) +
                         0.5 * rRight(i, j, k) * uRight(i, j, k).norm_squared();
      rELeft(i, j, k) = rLeft(i, j, k) * eLeft(i, j, k) +
                        0.5 * rLeft(i, j, k) * uLeft(i, j, k).norm_squared();
    }; // forall
  }
}

} // namespace hard::tasks::hydro

#endif // HARD_MODULE_HYDRO_RECONSTRUCT_HH
