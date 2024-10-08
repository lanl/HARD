#ifndef SPEC_UTILS_HH
#define SPEC_UTILS_HH

#include "mesh.hh"

#include <flecsi/flog.hh>

namespace spec::utils {

template<typename T>
FLECSI_INLINE_TARGET T
sqr(T t) {
  return t * t;
} // sq

template<std::size_t D>
inline mesh<D>::boundary_type
mesh_boundary(std::string const & b) {
  if(b == "flow") {
    return spec::mesh<D>::boundary_type::flow;
  }
  if(b == "reflecting") {
    return spec::mesh<D>::boundary_type::reflecting;
  }
  if(b == "periodic") {
    return spec::mesh<D>::boundary_type::periodic;
  }
  else {
    flog_fatal("invalid boundary type(" << b << ")");
  } // if
} // mesh_boundary

} // namespace spec::utils

#endif // SPEC_UTILS_HH
