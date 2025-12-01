#ifndef SPEC_LS_FIELDS_HH
#define SPEC_LS_FIELDS_HH

#include "../mesh.hh"
#include <flecsi/field/field.hpp>

namespace spec::level_set {
using flecsi::field;

// Signed distance-like scalar and integer material label per cell
// TODO: verify this
inline const field<double>::definition<spec::mesh<D>, is::cells> psi;
inline const field<int>::definition<spec::mesh<D>, is::cells> label;

inline const field<double>::definition<spec::mesh<D>, is::cells>
  Lwork; // HJ operator result
} // namespace spec::level_set
