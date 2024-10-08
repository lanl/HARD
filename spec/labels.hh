#ifndef SPEC_ENUM_HH
#define SPEC_ENUM_HH

namespace spec {

namespace ax {
/// Mesh axis.
/// The axis identifies the spatial axis in the Cartesian space.
enum axis {
  /// The x axis.
  x,
  /// The y axis.
  y,
  /// The z axis.
  z
};
} // namespace ax

namespace is {
enum index_space { cells };
}

namespace st {
enum dirs {
  /// The center stencil weight.
  c,
  /// The west stencil weight.
  w,
  /// The south stencil weight.
  s,
  /// The down stencil weight.
  d
};
}

namespace dm {
/// Mesh domains.
/// The domain identifies the supported iteration spaces on the mesh.
enum domain {
  /// This domain includes the locations of the unknowns of the problem state.
  quantities,
  /// This domain includes all mesh locations where it is possible to compute
  /// a 2nd-order slope approximation.
  predictor,
  /// This domain includes all mesh locations where it is possible to compute
  /// a flux across an interface from reconstructed predictor quantities.
  corrector,
  /// This domain includes all the degrees of freedom for the elliptic solver.
  interior,
  /// This domain includes all local cells, including halo and boundary cells.
  all,
  /// This domain includes all global cells.
  global
};
} // namespace dm

namespace bd {
/// Boundary.
/// Identifies the low or high boundary along a particular axis.
enum boundary {
  /// Low axis boundary.
  low,
  /// High axis boundary.
  high
}; // enum boundary

/// Boundary type.
/// The supported boundary types are derived from Leveque's
/// [Finite Volume Methods for Hyperbolc
/// Problems](https://www.amazon.com/Methods-Hyperbolic-Problems-Cambridge-Mathematics/dp/0521009243).
/// All boundary conditions are implemented using the \em ghost-cell approach
/// outlined in the text.
enum boundary_type {
  /// Zero-order extrapolation from the interior solution.
  flow,
  /// Cauchy problem solution for solid wall boundary.
  reflecting,
  /// Substitution with opposite interior solution.
  periodic
}; // enum boundary_type
} // namespace bd

} // namespace spec

#endif // SPEC_ENUM_HH
