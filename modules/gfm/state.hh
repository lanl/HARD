#ifndef HARD_MODULES_GFM_STATE_HH
#define HARD_MODULES_GFM_STATE_HH

namespace hard::gfm {
template<std::size_t D>
struct gfm {

  // Level set
  // Signed distance-like scalar and integer material label per cell
  // TODO: verify this
  static inline const field<double>::definition<mesh<D>, is::cells> psi;
  static inline const field<double>::definition<mesh<D>, is::cells> psi1, psi2,
    psi3;
  static inline const field<double>::definition<mesh<D>, is::cells> frozen_psi;

  // NEW: RHS storage for upwind operator (k / Lout)
  static inline const field<double>::definition<mesh<D>, is::cells> rhs;

  // Rename these for clarity
  static inline const field<int>::definition<mesh<D>, is::cells> material;
  static inline const field<int>::definition<mesh<D>, is::cells>
    frozen_material;

#if 0
  static inline const field<int>::definition<mesh<D>, is::cells> label;
  static inline const field<int>::definition<mesh<D>, is::cells> frozen_label;
#endif
};
} // namespace hard::gfm
#endif
