#ifndef HARD_MODULE_HYDRO_RECONSTRUCT_HH
#define HARD_MODULE_HYDRO_RECONSTRUCT_HH

#include "../modules/common/tasks/utils.hh"
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
// corresponding conserved variables on faces, and store them into `*right`
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
reconstruct_primitives_f(flecsi::exec::accelerator s,
  std::size_t reconstruction_axis,
  typename mesh<Dim>::template accessor<ro> m,
  // cell-centered primitive varibles
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    typename faces<Dim>::template accessor<wo, na>>> cons_faces_a,
  std::vector<std::tuple<typename field<vec<Dim>>::template accessor<ro, na>,
    typename faces_vec<Dim>::template accessor<wo, na>>>
    cons_faces_vec_a) noexcept {

  using common::tasks::utils::get_mdiota_policy;
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

template<std::size_t Dim, typename Limiter>
void
reconstruct_primitives(flecsi::exec::accelerator s,
  std::size_t reconstruction_axis,
  typename mesh<Dim>::template accessor<ro> m,
  // cell-centered primitive varibles
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    typename faces<Dim>::template accessor<wo, na>>> cons_faces_a,
  std::vector<std::tuple<typename field<vec<Dim>>::template accessor<ro, na>,
    typename faces_vec<Dim>::template accessor<wo, na>>>
    cons_faces_vec_a) noexcept {
  reconstruct_primitives_f<Dim, Limiter>(
    s, reconstruction_axis, m, cons_faces_a, cons_faces_vec_a);
}

template<std::size_t Dim, typename Limiter>
void
reconstruct_primitives_scalar(flecsi::exec::accelerator s,
  std::size_t reconstruction_axis,
  typename mesh<Dim>::template accessor<ro> m,
  // cell-centered primitive varibles
  std::vector<std::tuple<field<double>::accessor<ro, na>,
    typename faces<Dim>::template accessor<wo, na>>> cons_faces_a) noexcept {
  reconstruct_primitives_f<Dim, Limiter>(
    s, reconstruction_axis, m, cons_faces_a, {});
}

template<std::size_t Dim>
void
get_conserved(flecsi::exec::accelerator s,
  typename mesh<Dim>::template accessor<ro> m,
  typename faces<Dim>::template accessor<ro, na> rFace_a,
  typename faces_vec<Dim>::template accessor<ro, na> uFace_a,
  typename faces<Dim>::template accessor<ro, na> eFace_a,
  // Conserved variables on faces
  typename faces_vec<Dim>::template accessor<wo, na> ruFace_a,
  typename faces<Dim>::template accessor<wo, na> rEFace_a) noexcept {

  auto [r_right, r_left] = faces<Dim>::mdcolex(m, rFace_a);
  auto [u_right, u_left] = faces_vec<Dim>::mdcolex(m, uFace_a);
  auto [e_right, e_left] = faces<Dim>::mdcolex(m, eFace_a);
  auto [ru_right, ru_left] = faces_vec<Dim>::mdcolex(m, ruFace_a);
  auto [re_right, re_left] = faces<Dim>::mdcolex(m, rEFace_a);

  using common::tasks::utils::get_mdiota_policy;
  using spec::utils::sqr;

  if constexpr(Dim == 1) {

    s.executor().forall(i, (m.template cells<ax::x, dm::predictor>())) {

      // Compute conserved variables
      ru_right(i) = r_right(i) * u_right(i);
      ru_left(i) = r_left(i) * u_left(i);
      re_right(i) =
        r_right(i) * e_right(i) + 0.5 * r_right(i) * u_right(i).norm_squared();
      re_left(i) =
        r_left(i) * e_left(i) + 0.5 * r_left(i) * u_left(i).norm_squared();

    }; // forall
  }
  else if constexpr(Dim == 2) {

    auto mdpolicy_pp = get_mdiota_policy(r_right,
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(ji, mdpolicy_pp) {
      auto [j, i] = ji;

      // Compute conserved variables
      ru_right(i, j) = r_right(i, j) * u_right(i, j);
      ru_left(i, j) = r_left(i, j) * u_left(i, j);
      re_right(i, j) = r_right(i, j) * e_right(i, j) +
                       0.5 * r_right(i, j) * u_right(i, j).norm_squared();
      re_left(i, j) = r_left(i, j) * e_left(i, j) +
                      0.5 * r_left(i, j) * u_left(i, j).norm_squared();
    }; // forall
  }
  else { // Dim == 3

    auto mdpolicy_ppp = get_mdiota_policy(r_right,
      m.template cells<ax::z, dm::predictor>(),
      m.template cells<ax::y, dm::predictor>(),
      m.template cells<ax::x, dm::predictor>());

    s.executor().forall(kji, mdpolicy_ppp) {
      auto [k, j, i] = kji;

      // Compute conserved variables on faces
      ru_right(i, j, k) = r_right(i, j, k) * u_right(i, j, k);
      ru_left(i, j, k) = r_left(i, j, k) * u_left(i, j, k);
      re_right(i, j, k) =
        r_right(i, j, k) * e_right(i, j, k) +
        0.5 * r_right(i, j, k) * u_right(i, j, k).norm_squared();
      re_left(i, j, k) = r_left(i, j, k) * e_left(i, j, k) +
                         0.5 * r_left(i, j, k) * u_left(i, j, k).norm_squared();
    }; // forall
  }
}

} // namespace hard::tasks::hydro

#endif // HARD_MODULE_HYDRO_RECONSTRUCT_HH
