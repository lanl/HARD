#ifndef HARD_HYDRO_STATE_HH
#define HARD_HYDRO_STATE_HH

#include "hydro/state.hh"
#include "types.hh"

#include "hydro/tasks/utils.hh"

namespace hard {

/*----------------------------------------------------------------------------*
  Problem state.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct state : hydro::state<D> {

  flecsi::future<double> dtmin_;

  /*--------------------------------------------------------------------------*
    Topology slots.
    *--------------------------------------------------------------------------*/
  flecsi::topo::index::ptr ct; /* Color topology. */
  flecsi::topo::global::ptr gt;
  flecsi::topo::global::ptr dense_topology;

  // Deque for meshes in multigrid
  std::deque<typename mesh<D>::ptr> mh;

  // Define the fine grid
  typename mesh<D>::ptr & m = mh.emplace_back(typename mesh<D>::ptr());

  /*--------------------------------------------------------------------------*
    Color parameters (One per color using an index topology instance).
    *--------------------------------------------------------------------------*/

  /* Maximum characteristic speed for a color. */
  static inline const typename single<vec<D>>::template definition<index> lmax;

  static inline const single<double>::template definition<global> dt, t,
    dt_weighted;

  std::size_t lowest_level;
  std::size_t min_highest_level;
  std::size_t max_num_levels;

  spec::om::output_method output_method = spec::om::csv;

}; // struct state

} // namespace hard

#endif // HARD_STATE_HH
