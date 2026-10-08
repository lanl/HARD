#ifndef SPEC_TYPES_HH
#define SPEC_TYPES_HH

#include "labels.hh"
#include <flecsi/data.hh>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace tag {

using field_tags = unsigned;
enum field_tags_mask_t : field_tags {
  none = 0,
  boundary = 0x1 // receives boundary conditions
};

} // namespace tag

namespace spec {
inline constexpr flecsi::privilege na = flecsi::na, ro = flecsi::ro,
                                   wo = flecsi::wo, rw = flecsi::rw;

using flecsi::topo::global;
using flecsi::topo::index;
namespace py = pybind11;

using flecsi::field;
template<typename T>
using single = field<T, flecsi::data::single>;
using flecsi::data::multi;

// Wrapper arround FleCSI fields to attach a name to the fields
template<typename T, typename Topo>
struct tfield {
  using value_type = T;

  static constexpr std::size_t components = [] {
    if constexpr(std::is_same_v<T, double>)
      return std::size_t{1};
    else
      return T::Dim;
  }();

  explicit tfield(std::string n, tag::field_tags t = tag::none)
    : name(std::move(n)), tags(t) {}

  bool has(tag::field_tags t) const {
    return (tags & t) == t;
  }

  std::string name;
  tag::field_tags tags;
  typename field<T>::template definition<Topo, is::cells> def;

  template<typename... Args>
  auto operator()(Args &&... args) const
    -> decltype(def(std::forward<Args>(args)...)) {
    return def(std::forward<Args>(args)...);
  }
};

// Call f on every field of the tuple that carries all bits of tg.
template<typename Tuple, typename F>
void
for_each_field(Tuple const & t, tag::field_tags tg, F && f) {
  std::apply(
    [&](auto const &... fld) { ((fld.has(tg) ? (void)f(fld) : void()), ...); },
    t);
}

// ----------------------------------------------------------------------------
//  vec<D> type definitions
// ----------------------------------------------------------------------------
template<std::size_t D>
// requires(D <= 3)
struct vec {};

template<>
struct vec<1> {

  constexpr static std::size_t Dim = 1;

  FLECSI_INLINE_TARGET vec() = default;
  FLECSI_INLINE_TARGET ~vec() = default;

  FLECSI_INLINE_TARGET vec(double x) : vx(x) {}
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
  FLECSI_INLINE_TARGET double dot(const vec<1> & w) const {
    return vx * w[0];
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

  FLECSI_INLINE_TARGET vec(double x) {
    v[0] = x;
    v[1] = x;
  }
  FLECSI_INLINE_TARGET vec(double x, double y) {
    v[0] = x;
    v[1] = y;
  }

  FLECSI_INLINE_TARGET vec() = default;
  FLECSI_INLINE_TARGET ~vec() = default;

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
  FLECSI_INLINE_TARGET double dot(const vec<2> & w) const {
    return v[0] * w[0] + v[1] * w[1];
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

  FLECSI_INLINE_TARGET vec(double x) {
    v[0] = x;
    v[1] = x;
    v[2] = x;
  }
  FLECSI_INLINE_TARGET vec(double x, double y, double z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  FLECSI_INLINE_TARGET vec() = default;
  FLECSI_INLINE_TARGET ~vec() = default;

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

  FLECSI_INLINE_TARGET double dot(const vec<3> & w) const {
    return v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
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
FLECSI_INLINE_TARGET vec<D> &
operator*=(vec<D> & a, double s) {
  a[0] *= s;
  if constexpr(D > 1)
    a[1] *= s;
  if constexpr(D > 2)
    a[2] *= s;

  return a;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D> &
operator*=(vec<D> & a, const vec<D> & b) {
  a[0] *= b[0];
  if constexpr(D > 1)
    a[1] *= b[1];
  if constexpr(D > 2)
    a[2] *= b[2];

  return a;
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

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
min(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result[0] = std::min(a[0], b[0]);
  if constexpr(D > 1)
    result[1] = std::min(a[1], b[1]);
  if constexpr(D > 2)
    result[2] = std::min(a[2], b[2]);
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
max(const vec<D> & a, const vec<D> & b) {
  vec<D> result;
  result[0] = std::max(a[0], b[0]);
  if constexpr(D > 1)
    result[1] = std::max(a[1], b[1]);
  if constexpr(D > 2)
    result[2] = std::max(a[2], b[2]);
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
copysign(const vec<D> & mag, const vec<D> & sgn) {
  vec<D> result;
  result[0] = std::copysign(mag[0], sgn[0]);
  if constexpr(D > 1)
    result[1] = std::copysign(mag[1], sgn[1]);
  if constexpr(D > 2)
    result[2] = std::copysign(mag[2], sgn[2]);
  return result;
}

template<std::size_t D>
FLECSI_INLINE_TARGET vec<D>
copysign(double mag, const vec<D> & sgn) {
  vec<D> result;
  result[0] = std::copysign(mag, sgn[0]);
  if constexpr(D > 1)
    result[1] = std::copysign(mag, sgn[1]);
  if constexpr(D > 2)
    result[2] = std::copysign(mag, sgn[2]);
  return result;
}

// Component-wise minmax with scalar
// Returns pair of (min, max) for each component between scalar and vector
template<std::size_t D>
FLECSI_INLINE_TARGET std::tuple<vec<D>, vec<D>>
minmax(double scalar, const vec<D> & v) {
  vec<D> vmin, vmax;
  vmin[0] = std::min(scalar, v[0]);
  vmax[0] = std::max(scalar, v[0]);
  if constexpr(D > 1) {
    vmin[1] = std::min(scalar, v[1]);
    vmax[1] = std::max(scalar, v[1]);
  }
  if constexpr(D > 2) {
    vmin[2] = std::min(scalar, v[2]);
    vmax[2] = std::max(scalar, v[2]);
  }
  return {vmin, vmax};
}

// ----------------------------------------------------------------------------
//  vec<D> view type
// ----------------------------------------------------------------------------

template<std::size_t D>
struct vec_view {
  double * data;

  FLECSI_INLINE_TARGET double & operator[](std::size_t i) {
    return data[i];
  }

  FLECSI_INLINE_TARGET const double & operator[](std::size_t i) const {
    return data[i];
  }

  FLECSI_INLINE_TARGET vec_view & operator=(const vec<D> & v) {
    for(std::size_t d = 0; d < D; ++d)
      data[d] = v[d];
    return *this;
  }

  FLECSI_INLINE_TARGET operator vec<D>() const {
    vec<D> v{};
    for(std::size_t d = 0; d < D; ++d)
      v[d] = data[d];
    return v;
  }
};

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

struct __attribute__((visibility("default"))) config_py {
  config_py(const config_py &) = delete;
  config_py(config_py &&) = delete;
  config_py & operator=(const config_py &) = delete;
  config_py & operator=(config_py &&) = delete;

  config_py(const std::string & path) {
    py::gil_scoped_acquire _;
    py::object scope = py::globals();
    py::eval_file(path, scope);
    config_ = std::make_unique<py::object>(scope["config"]);
  }
  config_py(const py::object & object) {
    config_ = std::make_unique<py::object>(object);
  }

  bool contains(const std::string & key) const {
    py::gil_scoped_acquire _;
    return config_->contains(key);
  }

  // String key
  config_py operator[](const std::string & key) const {
    py::gil_scoped_acquire _;
    return config_py((*config_)[key.c_str()]);
  }

  // Index or integer key
  config_py operator[](const int & key) const {
    py::gil_scoped_acquire _;
    return config_py((*config_)[py::int_(key)]);
  }

  bool is_sequence() const {
    py::gil_scoped_acquire _;
    return py::isinstance<py::sequence>(*config_);
  }

  bool empty() const {
    return size() == 0;
  }

  template<typename T>
  T cast() const {
    py::gil_scoped_acquire _;
    return (*config_).cast<T>();
  }

  template<typename T>
  bool is() const {
    py::gil_scoped_acquire _;
    return py::isinstance<T>(*config_);
  }

  std::size_t size() const {
    py::gil_scoped_acquire _;
    return py::len(*config_);
  }

  // Destroy the unique_ptr in case of race condition
  ~config_py() {
    py::gil_scoped_acquire _;
    config_.reset();
  }

private:
  std::unique_ptr<py::object> config_;
}; // struct config_py
} // namespace spec

#endif // SPEC_TYPES_HH
