// spec/levelset/ls_fields.hh
#pragma once
#include <flecsi/field/field.hpp>
//TODO : is this correct to call?
#include "spec/mesh.hh" 

namespace spec::level_set {
  using flecsi::field;

  // Signed distance-like scalar and integer material label per cell
  // TODO: verify this
  inline const field<double>::definition<hard::mesh, hard::mesh::cells> psi;
  inline const field<int   >::definition<hard::mesh, hard::mesh::cells> label;

  inline const field<double>::definition<hard::mesh, hard::mesh::cells> Lwork; // HJ operator result
} // namespace spec::level_set

