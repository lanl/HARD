#ifndef HARD_MODULES_GFM_STATE_HH
#define HARD_MODULES_GFM_STATE_HH

#include "types.hh"

namespace hard::gfm {
template<std::size_t D>
struct state {
  struct level_set {
    // Level set
    // Signed distance-like scalar and integer material label per cell
    // TODO: verify this
    static inline const field<double>::definition<mesh<D>, is::cells> psi;
    static inline const field<double>::definition<mesh<D>, is::cells> psi1,
      psi2, psi3;
    static inline const field<double>::definition<mesh<D>, is::cells>
      frozen_psi;

    // NEW: RHS storage for upwind operator (k / Lout)
    static inline const field<double>::definition<mesh<D>, is::cells> rhs;

    // Rename these for clarity
    static inline const field<int>::definition<mesh<D>, is::cells> material;
    static inline const field<int>::definition<mesh<D>, is::cells>
      frozen_material;
  } ls;
  struct extrapolation {
    // Gradient component and unit normal
    static inline const typename field<vec<D>>::template definition<mesh<D>,
      is::cells>
      dp; // dp

    static inline const typename field<vec<D>>::template definition<mesh<D>,
      is::cells>
      normal; // dp

    static inline const field<double>::definition<mesh<D>, is::cells> q_scalar;
    static inline const field<double>::definition<mesh<D>, is::cells>
      q_scalar_known;
    static inline const field<double>::definition<mesh<D>, is::cells> q_scalar1;
    static inline const field<double>::definition<mesh<D>, is::cells> q_scalar2;
    static inline const field<double>::definition<mesh<D>, is::cells>
      q_new_scalar;
    // RHS for extrapolation
    // TODO : Maybe don't need this?
    static inline const field<double>::definition<mesh<D>, is::cells> q_rhs,
      q_rhs1;
  } extp;
}; // struct state
} // namespace hard::gfm
#endif
