#ifndef SPEC_RUNTIME_HH
#define SPEC_RUNTIME_HH

#include <flecsi/runtime.hh>
#include <utility>

namespace spec {
namespace detail {
template<template<template<std::size_t> typename S, std::size_t> typename C,
  template<std::size_t> typename S,
  template<std::size_t> typename F,
  auto CP,
  std::size_t... DD>
std::tuple<decltype([] {
  constexpr auto D = DD;
  return typename C<S, D>::template action<F<D>::action, CP>{};
}())...>
register_action(std::index_sequence<DD...>) {
  return {};
}
} // namespace detail

template<template<template<std::size_t> typename S, std::size_t> typename C,
  template<std::size_t> typename S,
  template<std::size_t> typename F,
  auto CP>
auto
register_action() {
  return detail::register_action<C, S, F, CP>(std::index_sequence<1, 2, 3>{});
}

namespace detail {
template<typename T>
struct is_tuple : std::false_type {};
template<typename... TT>
struct is_tuple<std::tuple<TT...>> : std::true_type {};

template<typename A, typename B>
concept same_size =
  (is_tuple<std::decay_t<A>>::value && is_tuple<std::decay_t<B>>::value &&
    (std::tuple_size_v<std::decay_t<A>> ==
      std::tuple_size_v<std::decay_t<B>>)) ||
  !(is_tuple<A>::value && is_tuple<B>::value);

template<typename T1, typename T2, typename F>
void visit_tuples(T1 && t1, T2 && t2, F && f)
  requires same_size<T1, T2>;

template<typename T1, typename T2, typename F, std::size_t... II>
void
visit_tuples_impl(T1 && t1, T2 && t2, F && f, std::index_sequence<II...>) {
  (...,
    visit_tuples(std::get<II>(std::forward<T1>(t1)),
      std::get<II>(std::forward<T2>(t2)),
      f));
}

template<typename T1, typename T2, typename F>
void
visit_tuples(T1 && t1, T2 && t2, F && f)
  requires same_size<T1, T2>
{
  if constexpr(is_tuple<std::decay_t<T1>>::value &&
               is_tuple<std::decay_t<T2>>::value) {
    visit_tuples_impl(std::forward<T1>(t1),
      std::forward<T2>(t2),
      std::forward<F>(f),
      std::make_index_sequence<std::tuple_size_v<std::decay_t<T1>>>{});
  }
  else {
    f(std::forward<T1>(t1), std::forward<T2>(t2));
  }
}
} // namespace detail

template<typename T1, typename T2>
auto
add_dependency(T1 && a, T2 && b) {
  detail::visit_tuples(
    a, b, [](auto & up, auto const & down) { up.add(down); });
  return true;
}

namespace detail {
template<template<template<std::size_t> typename S, std::size_t> typename C,
  template<std::size_t> typename S,
  std::size_t... DD,
  typename... AA>
auto
dispatch(flecsi::runtime & r,
  std::size_t d,
  std::index_sequence<DD...>,
  AA &&... aa) {
  bool ret = false;
  std::initializer_list<int>({(
    d == DD ? (ret = r.control<C<S, DD>>(std::forward<AA>(aa)...)), 0 : 0)...});
  return ret;
}
} // namespace detail

template<template<template<std::size_t> typename S, std::size_t> typename C,
  template<std::size_t> typename S,
  typename... AA>
auto
dispatch(flecsi::runtime & r, std::size_t d, AA &&... aa) {
  return detail::dispatch<C, S>(
    r, d, std::index_sequence<1, 2, 3>{}, std::forward<AA>(aa)...);
}
} // namespace spec

#endif // SPEC_RUNTIME_HH
