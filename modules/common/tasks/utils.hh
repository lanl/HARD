#ifndef HARD_MODULES_COMMON_TASKS_UTILS_HH
#define HARD_MODULES_COMMON_TASKS_UTILS_HH

#include <flecsi/data.hh>

namespace common::tasks::utils {

using namespace flecsi;

// Utility function to display variables associated with a name
template<std::size_t D>
void
display(std::vector<std::tuple<std::string, field<double>::accessor<ro, ro>>>
    v_f) noexcept {
  for(auto & vv : v_f) {
    auto [s, v] = vv;
    std::cout << s << ":" << '\n';
    for(int i = 0; i < v.span().size(); ++i) {
      std::cout << v[i] << '\n';
    }
  }
}

// Utility function to display min and max of a field
template<std::size_t D>
void
display_min_max(
  std::vector<std::tuple<std::string, field<double>::accessor<ro, ro>>>
    v_f) noexcept {
  for(auto & vv : v_f) {
    auto [s, v] = vv;
    auto [min, max] = std::ranges::minmax(v.span());
    std::cout << s << ": min,max = " << min << ", " << max << '\n';
  }
}

template<class M, typename IT>
FLECSI_INLINE_TARGET auto
get_mdiota_policy(const M & m, const IT & it1, const IT & it2) {
  using flecsi::exec::mdiota_view;
  using flecsi::exec::sub_range;

  auto b1 = *it1.begin();
  auto e1 = *it1.end();

  auto b2 = *it2.begin();
  auto e2 = *it2.end();

  return mdiota_view(m, sub_range{b1, e1}, sub_range{b2, e2});
} // get_mdiota_policy

template<class M, typename IT>
FLECSI_INLINE_TARGET auto
get_mdiota_policy(const M & m, const IT & it1, const IT & it2, const IT & it3) {
  using flecsi::exec::mdiota_view;
  using flecsi::exec::sub_range;

  auto b1 = *it1.begin();
  auto e1 = *it1.end();

  auto b2 = *it2.begin();
  auto e2 = *it2.end();

  auto b3 = *it3.begin();
  auto e3 = *it3.end();

  return mdiota_view(
    m, sub_range{b1, e1}, sub_range{b2, e2}, sub_range{b3, e3});
} // get_mdiota_policy

} // namespace common::tasks::utils

#endif // HARD_MODULES_COMMON_TASKS_UTILS_HH
