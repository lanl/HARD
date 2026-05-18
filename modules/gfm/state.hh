#ifndef HARD_MODULES_GFM_STATE_HH
#define HARD_MODULES_GFM_STATE_HH

#include "../modules/spec/eos.hh"
#include "tasks/utils.hh"
#include "types.hh"

namespace hard::gfm {

template<std::size_t D>
struct state {

  // Multi-mat EOS
  std::vector<eos::eos_wrapper> eos_mm;

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

    // Dynamics psi during extrapolation steps
    static inline const field<double>::definition<mesh<D>, is::cells> extp_psi;

    static inline const field<
      hard::tasks::util::extp_q<D>>::template definition<mesh<D>, is::cells>
      q_scalar; // scratch state being advanced by extrapolation PDE
    static inline const field<
      hard::tasks::util::extp_q<D>>::template definition<mesh<D>, is::cells>
      q_scalar_known; // fixed known-side values used to preserve psi<=0 region
    static inline const field<
      hard::tasks::util::extp_q<D>>::template definition<mesh<D>, is::cells>
      q_scalar1; // RK stage-1 scratch
    static inline const field<
      hard::tasks::util::extp_q<D>>::template definition<mesh<D>, is::cells>
      q_ext_plus; // RK final extrapolated state used by GFM on plus/right state
    static inline const field<
      hard::tasks::util::extp_q<D>>::template definition<mesh<D>, is::cells>
      q_ext_minus; // RK final extrapolated state used by GMF on minus/left
                   // state
    // RHS for extrapolation
    // TODO : Maybe don't need this?
    static inline const field<
      hard::tasks::util::extp_q<D>>::template definition<mesh<D>, is::cells>
      q_rhs, q_rhs1;

    // Normal faces for the surfaces
    struct faces_rl {
      faces_vec<D> normal_face;
    } f;

  } extp;
}; // struct state
} // namespace hard::gfm
#endif
