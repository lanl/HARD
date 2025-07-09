#ifndef SPEC_TYPES_HH
#define SPEC_TYPES_HH

#include "labels.hh"
#include <flecsi/data.hh>

namespace spec {
inline constexpr flecsi::privilege na = flecsi::na, ro = flecsi::ro,
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
// requires(D <= 3)
struct vec {};

template<>
struct vec<1> {

  constexpr static std::size_t Dim = 1;

  FLECSI_INLINE_TARGET vec<1>() = default;
  FLECSI_INLINE_TARGET ~vec<1>() = default;

  FLECSI_INLINE_TARGET vec<1>(double x) : vx(x) {}
  FLECSI_INLINE_TARGET double & operator[](std::size_t) {
    return vx;
  }
  FLECSI_INLINE_TARGET double operator[](std::size_t) const {
    return vx;
  }

  FLECSI_INLINE_TARGET double & x() {
    return vx;
  }
  FLECSI_INLINE_TARGET double x() const {
    return vx;
  }

  double vx;

  FLECSI_INLINE_TARGET double norm_squared() const {
    return vx * vx;
  }
  FLECSI_INLINE_TARGET double norm() const {
    return std::abs(vx);
  }
  FLECSI_INLINE_TARGET double & get([[maybe_unused]] const size_t idx) {
    assert(idx == 0 && "Invalid access for 1d vector");
    return vx;
  }
  FLECSI_INLINE_TARGET const double & get(
    [[maybe_unused]] const size_t idx) const {
    assert(idx == 0 && "Invalid access for 1d vector");

    return vx;
  }
};

template<>
struct vec<2> {

  constexpr static std::size_t Dim = 2;

  std::array<double, 2> v;

  FLECSI_INLINE_TARGET vec<2>(double x) {
    v[0] = x;
    v[1] = x;
  }
  FLECSI_INLINE_TARGET vec<2>(double x, double y) {
    v[0] = x;
    v[1] = y;
  }

  FLECSI_INLINE_TARGET vec<2>() = default;
  FLECSI_INLINE_TARGET ~vec<2>() = default;

  FLECSI_INLINE_TARGET double & operator[](std::size_t d) noexcept {
    return v[d];
  }

  FLECSI_INLINE_TARGET double operator[](std::size_t d) const noexcept {
    return v[d];
  }

  FLECSI_INLINE_TARGET double & x() {
    return v[0];
  }
  FLECSI_INLINE_TARGET double & y() {
    return v[1];
  }
  FLECSI_INLINE_TARGET double x() const {
    return v[0];
  }
  FLECSI_INLINE_TARGET double y() const {
    return v[1];
  }

  FLECSI_INLINE_TARGET double norm_squared() const {
    return v[0] * v[0] + v[1] * v[1];
  }
  FLECSI_INLINE_TARGET double norm() const {
    return std::sqrt(norm_squared());
  }
  FLECSI_INLINE_TARGET double & get(const size_t idx) {
    assert(idx <= 1 && "Invalid access for 2d vector");
    return v[idx];
  }
  FLECSI_INLINE_TARGET const double & get(const size_t idx) const {
    assert(idx <= 1 && "Invalid access for 2d vector");
    return v[idx];
  }
};

template<>
struct vec<3> {

  constexpr static std::size_t Dim = 3;

  std::array<double, 3> v;

  FLECSI_INLINE_TARGET vec<3>(double x) {
    v[0] = x;
    v[1] = x;
    v[2] = x;
  }
  FLECSI_INLINE_TARGET vec<3>(double x, double y, double z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  FLECSI_INLINE_TARGET vec<3>() = default;
  FLECSI_INLINE_TARGET ~vec<3>() = default;

  FLECSI_INLINE_TARGET double & operator[](std::size_t d) noexcept {
    return v[d];
  }

  FLECSI_INLINE_TARGET double operator[](std::size_t d) const noexcept {
    return v[d];
  }

  FLECSI_INLINE_TARGET double & x() {
    return v[0];
  }
  FLECSI_INLINE_TARGET double & y() {
    return v[1];
  }
  FLECSI_INLINE_TARGET double & z() {
    return v[2];
  }

  FLECSI_INLINE_TARGET double x() const {
    return v[0];
  }
  FLECSI_INLINE_TARGET double y() const {
    return v[1];
  }
  FLECSI_INLINE_TARGET double z() const {
    return v[2];
  }

  FLECSI_INLINE_TARGET double norm_squared() const {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }
  FLECSI_INLINE_TARGET double norm() const {
    return std::sqrt(norm_squared());
  }
  FLECSI_INLINE_TARGET double & get(const size_t idx) {
    assert(idx <= 2 && "Invalid access for 3d vector");
    return v[idx];
  }
  FLECSI_INLINE_TARGET const double & get(const size_t idx) const {
    assert(idx <= 2 && "Invalid access for 3d vector");
    return v[idx];
  }
};

// ----------------------------------------------------------------------------
//  Binary operators for vec<D> types
// ----------------------------------------------------------------------------

template<std::size_t D>
std::ostream &
operator<<(std::ostream & s, vec<D> const & v) {
  if constexpr(D == 1)
    s << "(" << v[0] << ")";
  else if constexpr(D == 2)
    s << "(" << v[0] << ", " << v[1] << ")";
  else
    s << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return s;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator+(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result[0] = a[0] + b[0];
  if constexpr(D > 1)
    result[1] = a[1] + b[1];
  if constexpr(D > 2)
    result[2] = a[2] + b[2];

  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator+(const vec<D> & a, double b) {
  vec<D> result;
  result[0] = a[0] + b;
  if constexpr(D > 1)
    result[1] = a[1] + b;
  if constexpr(D > 2)
    result[2] = a[2] + b;

  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator+(double b, const vec<D> & a) {
  vec<D> result;
  result[0] = a[0] + b;
  if constexpr(D > 1)
    result[1] = a[1] + b;
  if constexpr(D > 2)
    result[2] = a[2] + b;

  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator-(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result[0] = a[0] - b[0];
  if constexpr(D > 1)
    result[1] = a[1] - b[1];
  if constexpr(D > 2)
    result[2] = a[2] - b[2];
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET void
operator+=(vec<D> & a, const vec<D> & b) {
  a[0] += b[0];
  if constexpr(D > 1)
    a[1] += b[1];
  if constexpr(D > 2)
    a[2] += b[2];
}

template<std::size_t D>
FLECSI_INLINE_TARGET void
operator+=(vec<D> & a, double b) {
  a[0] += b;
  if constexpr(D > 1)
    a[1] += b;
  if constexpr(D > 2)
    a[2] += b;
}

template<std::size_t D>
FLECSI_INLINE_TARGET void
operator+=(double b, vec<D> & a) {
  a[0] += b;
  if constexpr(D > 1)
    a[1] += b;
  if constexpr(D > 2)
    a[2] += b;
}

template<std::size_t D>
FLECSI_INLINE_TARGET void
operator-=(vec<D> & a, const vec<D> & b) {
  a[0] -= b[0];
  if constexpr(D > 1)
    a[1] -= b[1];
  if constexpr(D > 2)
    a[2] -= b[2];
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator*(const vec<D> & a, double s) {
  vec<D> result;
  result[0] = a[0] * s;
  if constexpr(D > 1)
    result[1] = a[1] * s;
  if constexpr(D > 2)
    result[2] = a[2] * s;
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
  result[0] = a[0] * b[0];
  if constexpr(D > 1)
    result[1] = a[1] * b[1];
  if constexpr(D > 2)
    result[2] = a[2] * b[2];
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator/(const vec<D> & a, double s) {
  vec<D> result;
  result[0] = a[0] / s;
  if constexpr(D > 1)
    result[1] = a[1] / s;
  if constexpr(D > 2)
    result[2] = a[2] / s;
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
operator/(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result[0] = a[0] / b[0];
  if constexpr(D > 1)
    result[1] = a[1] / b[1];
  if constexpr(D > 2)
    result[2] = a[2] / b[2];
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
abs(const vec<D> & a) {
  vec<D> result;
  result[0] = std::abs(a[0]);
  if constexpr(D > 1)
    result[1] = std::abs(a[1]);
  if constexpr(D > 2)
    result[2] = std::abs(a[2]);
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
