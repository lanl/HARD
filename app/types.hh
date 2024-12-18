#ifndef HARD_TYPES_HH
#define HARD_TYPES_HH

#include <spec/control.hh>
#include <spec/types.hh>
#include <spec/utils.hh>

namespace hard {

/*----------------------------------------------------------------------------*
  Pull in some generic types.
 *----------------------------------------------------------------------------*/

template<template<std::size_t> typename S, std::size_t D>
using control = flecsi::run::control<spec::control_policy<S, D>>;

template<template<std::size_t> typename S, std::size_t D>
using control_policy = spec::control_policy<S, D>;

using cp = spec::cp;
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

inline constexpr flecsi::partition_privilege_t na = spec::na, ro = spec::ro,
                                               wo = spec::wo, rw = spec::rw;

/*----------------------------------------------------------------------------*
  Dual field.
 *----------------------------------------------------------------------------*/

template<typename T, std::size_t D>
struct dual_field {
  using type = const field<T>::template definition<mesh<D>, is::cells>;

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
  type fd_[2];
  int flip_{0};
};

} // namespace hard

#endif // HARD_TYPES_HH
