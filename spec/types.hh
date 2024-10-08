#ifndef SPEC_TYPES_HH
#define SPEC_TYPES_HH

#include "labels.hh"
#include <flecsi/data.hh>

namespace spec {
inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using flecsi::topo::global;
using flecsi::topo::index;

using flecsi::field;
template<typename T>
using single = field<T, flecsi::data::single>;
using flecsi::data::multi;

// ----------------------------------------------------------------------------
//  vec<D> type definitions
// ----------------------------------------------------------------------------
template<std::size_t D>
requires(D <= 3) struct vec {
};

template<>
struct vec<1> {
  union
  {
    struct {
      double x;
    };
    double components[1];
  };

  FLECSI_INLINE_TARGET double norm_squared() const {
    return x * x;
  }
  FLECSI_INLINE_TARGET double norm() const {
    return x;
  }
  FLECSI_INLINE_TARGET double & get(const size_t idx) {
#ifdef DEBUG
    if(idx != 0) {
      flog_fatal(
        "Invalid access occured vec<1>.get() function, trying to access "
        "the index "
        << idx);
    }
#endif
    return x;
  }
  FLECSI_INLINE_TARGET const double & get(const size_t idx) const {
#ifdef DEBUG
    if(idx != 0) {
      flog_fatal(
        "Invalid access occured vec<1>.get() function, trying to access "
        "the index "
        << idx);
    }
#endif
    return x;
  }
};

template<>
struct vec<2> {
  union
  {
    struct {
      double x;
      double y;
    };
    double components[2];
  };

  FLECSI_INLINE_TARGET double norm_squared() const {
    return x * x + y * y;
  }
  FLECSI_INLINE_TARGET double norm() const {
    return std::sqrt(norm_squared());
  }
  FLECSI_INLINE_TARGET double & get(const size_t idx) {
#ifdef DEBUG
    if((idx < 0) or (idx > 1)) {
      flog_fatal(
        "Invalid access occured vec<2>.get() function, trying to access "
        "the index "
        << idx);
    }
#endif
    return components[idx];
  }
  FLECSI_INLINE_TARGET const double & get(const size_t idx) const {
#ifdef DEBUG
    if((idx < 0) or (idx > 1)) {
      flog_fatal(
        "Invalid access occured vec<2>.get() function, trying to access "
        "the index "
        << idx);
    }
#endif
    return components[idx];
  }
};

template<>
struct vec<3> {
  union
  {
    struct {
      double x;
      double y;
      double z;
    };
    double components[3];
  };

  FLECSI_INLINE_TARGET double norm_squared() const {
    return x * x + y * y + z * z;
  }
  FLECSI_INLINE_TARGET double norm() const {
    return std::sqrt(norm_squared());
  }
  FLECSI_INLINE_TARGET double & get(const size_t idx) {
#ifdef DEBUG
    if((idx < 0) or (idx > 2)) {
      flog_fatal(
        "Invalid access occured vec<3>.get() function, trying to access "
        "the index "
        << idx);
    }
#endif
    return components[idx];
  }
  FLECSI_INLINE_TARGET const double & get(const size_t idx) const {
#ifdef DEBUG
    if((idx < 0) or (idx > 2)) {
      flog_fatal(
        "Invalid access occured vec<3>.get() function, trying to access "
        "the index "
        << idx);
    }
#endif
    return components[idx];
  }
};

// ----------------------------------------------------------------------------
//  Binary operators for vec<D> types
// ----------------------------------------------------------------------------

template<std::size_t D>
std::ostream &
operator<<(std::ostream & s, vec<D> const & v) {
  if constexpr(D == 1)
    s << "(" << v.x << ")";
  else if constexpr(D == 2)
    s << "(" << v.x << ", " << v.y << ")";
  else
    s << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return s;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator+(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result.x = a.x + b.x;
  if constexpr(D > 1)
    result.y = a.y + b.y;
  if constexpr(D > 2)
    result.z = a.z + b.z;

  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator-(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result.x = a.x - b.x;
  if constexpr(D > 1)
    result.y = a.y - b.y;
  if constexpr(D > 2)
    result.z = a.z - b.z;
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET void
operator+=(vec<D> & a, const vec<D> & b) {
  a.x += b.x;
  if constexpr(D > 1)
    a.y += b.y;
  if constexpr(D > 2)
    a.z += b.z;
}

template<std::size_t D>
FLECSI_INLINE_TARGET void
operator-=(vec<D> & a, const vec<D> & b) {
  a.x -= b.x;
  if constexpr(D > 1)
    a.y -= b.y;
  if constexpr(D > 2)
    a.z -= b.z;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator*(const vec<D> & a, double s) {
  vec<D> result;
  result.x = a.x * s;
  if constexpr(D > 1)
    result.y = a.y * s;
  if constexpr(D > 2)
    result.z = a.z * s;
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator*(double s, const vec<D> & a) {
  return a * s;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator*(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result.x = a.x * b.x;
  if constexpr(D > 1)
    result.y = a.y * b.y;
  if constexpr(D > 2)
    result.z = a.z * b.z;
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator/(const vec<D> & a, double s) {
  vec<D> result;
  result.x = a.x / s;
  if constexpr(D > 1)
    result.y = a.y / s;
  if constexpr(D > 2)
    result.z = a.z / s;
  return result;
}

// ----------------------------------------------------------------------------
//  Tensor types definitions
// ----------------------------------------------------------------------------

//
// @ Yoonsoo Kim (Jun 25) :
//  We can possibly make the following tensor struct more generic, for example,
//  by adding a method like `.get(i,j,k..)` returning the contained values.
//  However, rank 2 tensor is what we need so far (for radiation pressure), so
//  it would be enough for now to just write out all cases by hand..
//

// Declare enum for the tensor rank to avoid confusion with the dimension
// template parameter, since both of them are unsigned integers.
enum tensor_rank {
  Two,
  Three
  // ...
};

template<std::size_t Dim, tensor_rank Rank>
struct tensor {};

template<>
struct tensor<1, tensor_rank::Two> {
  double xx;
  friend std::ostream & operator<<(std::ostream & s, tensor const & t) {
    s << "(" << t.xx << ")";
    return s;
  }
};

template<>
struct tensor<2, tensor_rank::Two> {
  double xx, xy, yx, yy;
  friend std::ostream & operator<<(std::ostream & s, tensor const & t) {
    s << "(" << t.xx << ", " << t.xy << ", " << t.yx << ", " << t.yy << ")";
    return s;
  }
};

template<>
struct tensor<3, tensor_rank::Two> {
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;
  friend std::ostream & operator<<(std::ostream & s, tensor const & t) {
    s << "(" << t.xx << ", " << t.xy << ", " << t.xz << ", " << t.yx << ", "
      << t.yy << ", " << t.yz << ", " << t.zx << ", " << t.zy << ", " << t.zz
      << ")";
    return s;
  }
};

// ----------------------------------------------------------------------------

template<std::size_t D>
struct stencil {
  std::array<double, D + 1> weights_;
  FLECSI_INLINE_TARGET double operator[](st::dirs d) const {
    return weights_[d];
  }
  FLECSI_INLINE_TARGET double & operator[](st::dirs d) {
    return weights_[d];
  }
  friend std::ostream & operator<<(std::ostream & s, stencil const & ws) {
    s << "c: " << ws[st::c] << " w: " << ws[st::w];
    if constexpr(D == 2 || D == 3) {
      s << " s: " << ws[st::s];
    }
    if constexpr(D == 3) {
      s << " d: " << ws[st::d];
    }
    return s;
  }
};
} // namespace spec

#endif // SPEC_TYPES_HH
