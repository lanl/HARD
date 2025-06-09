#ifndef HARD_UTILS_HH
#define HARD_UTILS_HH

#include "types.hh"

namespace hard::util {

/*!
  Return the total number of colors across all axes from the given @em colors
  reference.

  @tparam D The spatial dimension.

  @param cd The @em colors data structure of axis colors.
 */
template<std::size_t D>
std::size_t
axes_colors(color_distribution const & cd) {
  if(D == 1) {
    return cd[0];
  }
  else if(D == 2) {
    return cd[0] * cd[1];
  }
  else /* D == 3 */ {
    return cd[0] * cd[1] * cd[2];
  }
} // colors

} // namespace hard::util

#endif // HARD_UTILS_HH
