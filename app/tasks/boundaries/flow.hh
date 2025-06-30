#ifndef TASKS_BOUNDARIES_FLOW_HH
#define TASKS_BOUNDARIES_FLOW_HH

#include "../../tasks/utils.hh"
#include "../../types.hh"

namespace hard::tasks {

using hard::tasks::util::bl;

template<int D>
struct flow {};

template<>
struct flow<1> {
  template<typename T>
  void operator()(flecsi::util::mdcolex<T, 1> a, int i, int level) const {
    if(level == bl::low)
      for(int m = 0; m < ghost_zone_size; ++m)
        a(m) = a(ghost_zone_size);

    if(level == bl::high)
      for(int m = 0; m < ghost_zone_size; ++m)
        a(i - 1 - m) = a(i - 1 - ghost_zone_size);
  }
  int ghost_zone_size;
};

template<>
struct flow<2> {
  template<typename T>
  void operator()(int axis,
    flecsi::util::mdcolex<T, 2> a,
    int i,
    int j,
    int level) const {
    if(level == bl::low) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(m, j) = a(ghost_zone_size, j);
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, m) = a(i, ghost_zone_size);
    }
    if(level == bl::high) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i - 1 - m, j) = a(i - 1 - ghost_zone_size, j);
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j - 1 - m) = a(i, j - 1 - ghost_zone_size);
    }
  }
  int ghost_zone_size;
};

template<>
struct flow<3> {
  template<typename T>
  void operator()(int axis,
    flecsi::util::mdcolex<T, 3> a,
    int i,
    int j,
    int k,
    int level) const {
    if(level == bl::low) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(m, j, k) = a(ghost_zone_size, j, k);
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, m, k) = a(i, ghost_zone_size, k);
      if(axis == ax::z)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j, m) = a(i, j, ghost_zone_size);
    }
    if(level == bl::high) {
      if(axis == ax::x)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i - 1 - m, j, k) = a(i - 1 - ghost_zone_size, j, k);
      if(axis == ax::y)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j - 1 - m, k) = a(i, j - 1 - ghost_zone_size, k);
      if(axis == ax::z)
        for(int m = 0; m < ghost_zone_size; ++m)
          a(i, j, k - 1 - m) = a(i, j, k - 1 - ghost_zone_size);
    }
  }
  int ghost_zone_size;
};
} // namespace hard::tasks
#endif
