#ifndef HARD_HYDRO_TYPES_HH
#define HARD_HYDRO_TYPES_HH

#include <hydro/types.hh>

#include <spec/control.hh>
#include <spec/eos.hh>
#include <spec/types.hh>
#include <spec/utils.hh>

#include <flecsi/utilities.hh>

#include "spec/vector_field.hh"

namespace hard {

/*----------------------------------------------------------------------------*
  Pull in some generic types.
 *----------------------------------------------------------------------------*/

template<template<std::size_t> typename S, std::size_t D>
using control = flecsi::run::control<spec::mmhydro::control_policy<S, D>>;

template<template<std::size_t> typename S, std::size_t D>
using control_policy = spec::mmhydro::control_policy<S, D>;
using color_distribution = flecsi::topo::narray_impl::colors;

using cp = spec::mmhydro::cp;
template<std::size_t D>
using mesh = spec::mesh<D>;
using spec::field;
using spec::global;
using spec::index;
using spec::multi;
using spec::single;
using spec::stencil;
using spec::vec;
using spec::st::dirs;

struct annotation_execution_action
  : flecsi::util::annotation::context<annotation_execution_action> {
  static constexpr char name[] = "HARD-Action";
};
using annotation_action_region =
  flecsi::util::annotation::region<annotation_execution_action>;

struct annotation_execution_task
  : flecsi::util::annotation::context<annotation_execution_task> {
  static constexpr char name[] = "HARD-Task";
};
using annotation_task_region =
  flecsi::util::annotation::region<annotation_execution_task>;

template<class Region>
using annotation_rguard = flecsi::util::annotation::rguard<Region>;

/*----------------------------------------------------------------------------*
  Namespace labels.
 *----------------------------------------------------------------------------*/

namespace utils = spec::utils;
namespace is = spec::is;
namespace ax = spec::ax;
namespace dm = spec::dm;
namespace bd = spec::bd;

// Suppress irritating warnings in clangd. Probably a better way to do this...
#define NS_WARN_SUPPRESS(n, m)                                                 \
  inline void n##_suppress() {                                                 \
    using m;                                                                   \
  }

NS_WARN_SUPPRESS(utils, utils::sqr);
NS_WARN_SUPPRESS(is, is::cells);
NS_WARN_SUPPRESS(ax, ax::x);
NS_WARN_SUPPRESS(dm, dm::all);
NS_WARN_SUPPRESS(bd, bd::low);

/*----------------------------------------------------------------------------*
  More concise privileges.
 *----------------------------------------------------------------------------*/

inline constexpr flecsi::privilege na = spec::na, ro = spec::ro, wo = spec::wo,
                                   rw = spec::rw;

/*----------------------------------------------------------------------------*
  Dual field.
 *----------------------------------------------------------------------------*/

template<typename T, std::size_t D>
struct dual_field {
  using type = const typename field<T>::template definition<mesh<D>, is::cells>;

  auto flip() {
    ++flip_;
  } // flip

  /*!
    Return the requested field reference.

    @tparam S The topology slot type.
    @param  s The topology slot instance.
    @param  i The index of the field value to return.
   */
  template<typename S>
  auto operator()(S & s, int i = 0) const {
    return fd_[(flip_ + i) % 2](s);
  } // operator

private:
  type fd_[2] = {};
  int flip_{0};
};

template<std::size_t D>
struct faces {

  using ref = flecsi::field<double>::Reference<mesh<D>, is::cells>;
  using ref_pair = std::tuple<ref, ref>;

  std::tuple<field<double>::definition<mesh<D>, is::cells>, // right
    field<double>::definition<mesh<D>, is::cells> // left
    >
    f;

  template<flecsi::privilege P1, flecsi::privilege P2>
  using accessor = std::tuple<field<double>::accessor<P1, P2>,
    field<double>::accessor<P1, P2>>;

  auto operator()(const mesh<D>::ptr & s) {
    return std::make_tuple(std::get<0>(f)(*s), std::get<1>(f)(*s));
  }

  template<class MeshAcc, class Tuple>
  static auto mdcolex(MeshAcc && m, Tuple & a) {
    return std::apply(
      [&](auto &... x) {
        return std::tuple{m.template mdcolex<is::cells>(x)...};
      },
      a);
  }
};

template<std::size_t D>
struct faces_vec {

  std::tuple<vector_field_double<D>, // right
    vector_field_double<D> // left
    >
    f;

  template<flecsi::privilege P1, flecsi::privilege P2>
  using accessor =
    std::tuple<typename vector_field_double<D>::template accessor<P1, P2>,
      typename vector_field_double<D>::template accessor<P1, P2>>;

  auto operator()(const mesh<D>::ptr & s) {
    return std::make_tuple(std::get<0>(f)(*s), std::get<1>(f)(*s));
  }

  template<class MeshAcc, class Tuple>
  static auto mdcolex(MeshAcc && m, Tuple & a) {
    return std::apply(
      [&](auto & vec_refs_right, auto & vec_refs_left) {
        auto right_accessor =
          vector_field_double<D>::mdcolex(m, vec_refs_right);
        auto left_accessor = vector_field_double<D>::mdcolex(m, vec_refs_left);

        return std::make_tuple(
          std::move(right_accessor), std::move(left_accessor));
      },
      a);
  }
};

template<std::size_t D>
struct RK {
  std::tuple<field<double>::definition<mesh<D>, is::cells>, // mass_density
    field<double>::definition<mesh<D>, is::cells>, // total_energy_density
    vector_field_double<D> // momentum_energy_density
    >
    f;

  template<flecsi::privilege P1, flecsi::privilege P2>
  using accessor = std::tuple<field<double>::accessor<P1, P2>,
    field<double>::accessor<P1, P2>,
    typename vector_field_double<D>::template accessor<P1, P2>>;

  auto & mass_density() {
    return std::get<0>(f);
  }

  const auto & mass_density() const {
    return std::get<0>(f);
  }

  auto & total_energy_density() {
    return std::get<1>(f);
  }

  const auto & total_energy_density() const {
    return std::get<1>(f);
  }

  auto & momentum_energy_density() {
    return std::get<2>(f);
  }

  const auto & momentum_energy_density() const {
    return std::get<2>(f);
  }

  auto operator()(const mesh<D>::ptr & s) {
    return std::make_tuple(
      std::get<0>(f)(*s), std::get<1>(f)(*s), std::get<2>(f)(*s));
  }

  template<class MeshAcc, class Tuple>
  static auto mdcolex(MeshAcc && m, Tuple & a) {
    return std::apply(
      [&](auto & mass_ref, auto & energy_ref, auto & momentum_ref) {
        auto mass_acc = m.template mdcolex<is::cells>(mass_ref);
        auto energy_acc = m.template mdcolex<is::cells>(energy_ref);
        auto momentum_acc = vector_field_double<D>::mdcolex(m, momentum_ref);

        return std::make_tuple(
          std::move(mass_acc), std::move(energy_acc), std::move(momentum_acc));
      },
      a);
  }
};

} // namespace hard

#endif // HARD_HYDRO_TYPES_HH
