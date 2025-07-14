#ifndef TASKS_BOUNDARIES_REFLECTIVE_HH
#define TASKS_BOUNDARIES_REFLECTIVE_HH

#include "../../tasks/utils.hh"
#include "../../types.hh"

namespace hard::tasks {

using hard::tasks::util::bl;

template<size_t D>
struct reflective {};

template<>
struct reflective<1> {

  FLECSI_INLINE_TARGET reflective<1>(int gzs) : ghost_zone_size(gzs){};

  template<typename T>
  FLECSI_INLINE_TARGET void
  operator()(flecsi::util::mdcolex<T, 1> a, int i, int level) const {
    if(level == bl::low) {
      for(int m = 0; m < ghost_zone_size; ++m)
        if constexpr(std::is_same_v<T, double>)
          a(m) = a(2 * ghost_zone_size - 1 - m);
        else
          a(m) = -1 * a(2 * ghost_zone_size - 1 - m);
    }
    if(level == bl::high) {
      for(int m = 0; m < ghost_zone_size; ++m)
        if constexpr(std::is_same_v<T, double>)
          a(i - 1 - m) = a(i - 2 * ghost_zone_size + m);
        else
          a(i - 1 - m) = -1 * a(i - 2 * ghost_zone_size + m);
    }
  }
  int ghost_zone_size;
};

template<>
struct reflective<2> {

  FLECSI_INLINE_TARGET reflective<2>(int gzs) : ghost_zone_size(gzs){};

  template<typename T>
  FLECSI_INLINE_TARGET void operator()(int axis,
    flecsi::util::mdcolex<T, 2> a,
    int i,
    int j,
    int level) const {
    if(level == bl::low) {
      if(axis == ax::x) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(m, j) = a(2 * ghost_zone_size - 1 - m, j);
          else {
            a(m, j).x() = -1 * a(2 * ghost_zone_size - 1 - m, j).x();
            a(m, j).y() = a(2 * ghost_zone_size - 1 - m, j).y();
          }
      }
      if(axis == ax::y) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i, m) = a(i, 2 * ghost_zone_size - 1 - m);
          else {
            a(i, m).x() = a(i, 2 * ghost_zone_size - 1 - m).x();
            a(i, m).y() = -1 * a(i, 2 * ghost_zone_size - 1 - m).y();
          }
      }
    }
    if(level == bl::high) {
      if(axis == ax::x) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i - 1 - m, j) = a(i - 2 * ghost_zone_size + m, j);
          else {
            a(i - 1 - m, j).x() = -1 * a(i - 2 * ghost_zone_size + m, j).x();
            a(i - 1 - m, j).y() = a(i - 2 * ghost_zone_size + m, j).y();
          }
        if(axis == ax::y) {
          for(int m = 0; m < ghost_zone_size; ++m)
            if constexpr(std::is_same_v<T, double>)
              a(i, j - 1 - m) = a(i, j - 2 * ghost_zone_size + m);
            else {
              a(i, j - 1 - m).x() = a(i, j - 2 * ghost_zone_size + m).y();
              a(i, j - 1 - m).y() = -1 * a(i, j - 2 * ghost_zone_size + m).y();
            }
        }
      }
    }
  }
  int ghost_zone_size;
};

template<>
struct reflective<3> {

  FLECSI_INLINE_TARGET reflective<3>(int gzs) : ghost_zone_size(gzs){};

  template<typename T>
  FLECSI_INLINE_TARGET void operator()(int axis,
    flecsi::util::mdcolex<T, 3> a,
    int i,
    int j,
    int k,
    int level) const {
    if(level == bl::low) {
      if(axis == ax::x) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(m, j, k) = a(2 * ghost_zone_size - 1 - m, j, k);
          else {
            a(m, j, k).x() = -1 * a(2 * ghost_zone_size - 1 - m, j, k).x();
            a(m, j, k).y() = a(2 * ghost_zone_size - 1 - m, j, k).y();
            a(m, j, k).z() = a(2 * ghost_zone_size - 1 - m, j, k).z();
          }
      }
      if(axis == ax::y) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i, m, k) = a(i, 2 * ghost_zone_size - 1 - m, k);
          else {
            a(i, m, k).x() = a(i, 2 * ghost_zone_size - 1 - m, k).x();
            a(i, m, k).y() = -1 * a(i, 2 * ghost_zone_size - 1 - m, k).y();
            a(i, m, k).z() = a(i, 2 * ghost_zone_size - 1 - m, k).z();
          }
      }
      if(axis == ax::z) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i, j, m) = a(i, j, 2 * ghost_zone_size - 1 - m);
          else {
            a(i, j, m).x() = a(i, j, 2 * ghost_zone_size - 1 - m).x();
            a(i, j, m).y() = a(i, j, 2 * ghost_zone_size - 1 - m).y();
            a(i, j, m).z() = -1 * a(i, j, 2 * ghost_zone_size - 1 - m).z();
          }
      }
    }
    if(level == bl::high) {
      if(axis == ax::x) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i - 1 - m, j, k) = a(i - 2 * ghost_zone_size + m, j, k);
          else {
            a(i - 1 - m, j, k).x() =
              -1 * a(i - 2 * ghost_zone_size + m, j, k).x();
            a(i - 1 - m, j, k).y() = a(i - 2 * ghost_zone_size + m, j, k).y();
            a(i - 1 - m, j, k).z() = a(i - 2 * ghost_zone_size + m, j, k).z();
          }
      }
      if(axis == ax::y) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i, j - 1 - m, k) = a(i, j - 2 * ghost_zone_size + m, k);
          else {
            a(i, j - 1 - m, k).x() = a(i, j - 2 * ghost_zone_size + m, k).x();
            a(i, j - 1 - m, k).y() =
              -1 * a(i, j - 2 * ghost_zone_size + m, k).y();
            a(i, j - 1 - m, k).z() = a(i, j - 2 * ghost_zone_size + m, k).z();
          }
      }
      if(axis == ax::z) {
        for(int m = 0; m < ghost_zone_size; ++m)
          if constexpr(std::is_same_v<T, double>)
            a(i, j, k - 1 - m) = a(i, j, k - 2 * ghost_zone_size + m);
          else {
            a(i, j, k - 1 - m).x() = a(i, j, k - 2 * ghost_zone_size + m).x();
            a(i, j, k - 1 - m).y() = a(i, j, k - 2 * ghost_zone_size + m).y();
            a(i, j, k - 1 - m).z() =
              -1 * a(i, j, k - 2 * ghost_zone_size + m).z();
          }
      }
    }
  }
  int ghost_zone_size;
};
} // namespace hard::tasks
#endif
