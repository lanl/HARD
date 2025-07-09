#ifndef TASKS_BOUNDARIES_DIRICHLET_HH
#define TASKS_BOUNDARIES_DIRICHLET_HH

#include "../../tasks/utils.hh"
#include "../../types.hh"

namespace hard::tasks {

using hard::tasks::util::bl;

template<size_t D>
struct dirichlet {};

template<>
struct dirichlet<1> {

  FLECSI_INLINE_TARGET dirichlet<1>(int gzs) : ghost_zone_size(gzs) {}

  template<typename T>
  FLECSI_INLINE_TARGET void
  operator()(flecsi::util::mdcolex<T, 1> a, double v, int i, int level) const {
    if(level == bl::low)
      for(int m = 0; m < ghost_zone_size; ++m)
        a(m) = v;
    if(level == bl::high)
      for(int m = 0; m < ghost_zone_size; ++m)
        a(i - 1 - m) = v;
  }
  int ghost_zone_size;
};

template<>
struct dirichlet<2> {

  FLECSI_INLINE_TARGET dirichlet<2>(int gzs) : ghost_zone_size(gzs) {}

  template<typename T>
  FLECSI_INLINE_TARGET void operator()(int axis,
    flecsi::util::mdcolex<T, 2> a,
    double v,
    int i,
    int j,
    int level) const {
    if(level == bl::low) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(m, j) = v;
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, m) = v;
    }
    if(level == bl::high) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i - 1 - m, j) = v;
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j - 1 - m) = v;
    }
  }
  int ghost_zone_size;
};

template<>
struct dirichlet<3> {

  FLECSI_INLINE_TARGET dirichlet<3>(int gzs) : ghost_zone_size(gzs) {}

  template<typename T>
  FLECSI_INLINE_TARGET void operator()(int axis,
    flecsi::util::mdcolex<T, 3> a,
    double v,
    int i,
    int j,
    int k,
    int level) const {
    if(level == bl::low) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(m, j, k) = v;
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, m, k) = v;
      if(axis == ax::z)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j, m) = v;
    }
    if(level == bl::high) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i - 1 - m, j, k) = v;
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j - 1 - m, k) = v;
      if(axis == ax::z)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j, k - 1 - m) = v;
    }
  }
  int ghost_zone_size;
};
} // namespace hard::tasks
#endif
