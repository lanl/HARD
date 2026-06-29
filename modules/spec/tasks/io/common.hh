#ifndef HARD_TASKS_IO_COMMON_HH
#define HARD_TASKS_IO_COMMON_HH

#include "state.hh"
#include <algorithm>
#include <string>
#include <vector>

namespace hard::tasks::io {

// Helper to extract field accessors and names from tuple vectors
template<std::size_t D>
struct field_extractor {
  std::vector<flecsi::util::mdcolex<const double, D>> fields;
  std::vector<flecsi::util::mdcolex<const vec<D>, D>> fields_vectors;
  std::vector<std::string> field_name;
  std::vector<std::string> field_vector_name;

  field_extractor(const auto & m,
    uint32_t i,
    const std::vector<std::tuple<multi<field<double>::accessor<ro, ro>>,
      std::string>> & field_ma,
    const std::vector<
      std::tuple<multi<typename field<vec<D>>::template accessor<ro, ro>>,
        std::string>> & field_vector_ma) {

    std::transform(field_ma.begin(),
      field_ma.end(),
      std::back_inserter(fields),
      [&m, &i](const auto & tuple) -> flecsi::util::mdcolex<const double, D> {
        return m.template mdcolex<is::cells>(std::get<0>(tuple).accessors()[i]);
      });
    std::transform(field_ma.begin(),
      field_ma.end(),
      std::back_inserter(field_name),
      [](const auto & tuple) { return std::get<1>(tuple); });
    std::transform(field_vector_ma.begin(),
      field_vector_ma.end(),
      std::back_inserter(fields_vectors),
      [&m, &i](const auto & tuple) -> flecsi::util::mdcolex<const vec<D>, D> {
        return m.template mdcolex<is::cells>(std::get<0>(tuple).accessors()[i]);
      });
    std::transform(field_vector_ma.begin(),
      field_vector_ma.end(),
      std::back_inserter(field_vector_name),
      [](const auto & tuple) { return std::get<1>(tuple); });
  }
};

} // namespace hard::tasks::io

#endif // HARD_TASKS_IO_COMMON_HH
