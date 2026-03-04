#ifndef SPEC_RUNTIME_HH
#define SPEC_RUNTIME_HH

#include <flecsi/runtime.hh>

// clang-format off
namespace detail {
template<
  template<
    template<std::size_t> typename S,
    std::size_t
  > typename C,
  template<std::size_t> typename S,
  std::size_t... DD,
  typename... AA
>
auto
dispatch(flecsi::runtime & r,
  std::size_t d,
  std::index_sequence<DD...>,
  AA &&... aa) {
  bool ret = false;
  std::initializer_list<int>({(
    d == DD ? (ret = r.control<C<S, DD>>(std::forward<AA>(aa)...)), 0 : 0)...});
  return ret;
} // dispacth
} // namespace detail

template<
  template<
    template<std::size_t> typename S,
    std::size_t
  > typename C,
  template<std::size_t>
  typename S,
  typename... AA>
auto
dispatch(flecsi::runtime & r, std::size_t d, AA &&... aa) {
  return detail::dispatch<C, S>(
    r, d, std::index_sequence<1, 2, 3>{}, std::forward<AA>(aa)...);
} // dispacth
// clang-format on

#endif // SPEC_RUNTIME_HH
