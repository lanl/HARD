#ifndef SPEC_RECONSTRUCT_HH
#define SPEC_RECONSTRUCT_HH

#include "types.hh"
#include "utils.hh"

namespace spec::reconstructors {

namespace detail {

// size 3 gives {-1, 0, 1} and size 5 gives {-2, -1, 0, 1, 2}.
// This runs at compile time
template<std::size_t N>
FLECSI_INLINE_TARGET constexpr std::array<int, N>
centered_offsets() {
  std::array<int, N> offs{};
  for(std::size_t i = 0; i < N; ++i)
    offs[i] = static_cast<int>(i) - static_cast<int>(N / 2);
  return offs;
}

// CRTP base implementing the axis/index dispatch shared by every
// reconstructor.
template<typename Derived, std::size_t N>
struct stencil {

  static constexpr std::array<int, N> offsets = centered_offsets<N>();

  // 1D
  template<typename A, typename T>
  FLECSI_INLINE_TARGET static auto reconstruct(const T & i, const A & acc) {
    return gather(acc, i, std::make_index_sequence<N>{});
  }

  // 2D
  template<typename A, typename T>
  FLECSI_INLINE_TARGET static auto
  reconstruct(const int & axis, const T & i, const T & j, const A & acc) {
    return gather(acc, axis, i, j, std::make_index_sequence<N>{});
  }

  // 3D
  template<typename A, typename T>
  FLECSI_INLINE_TARGET static auto reconstruct(const int & axis,
    const T & i,
    const T & j,
    const T & k,
    const A & acc) {
    return gather(acc, axis, i, j, k, std::make_index_sequence<N>{});
  }

private:
  // 1D
  template<typename A, typename T, std::size_t... K>
  FLECSI_INLINE_TARGET static auto
  gather(const A & acc, const T & i, std::index_sequence<K...>) {
    return Derived::stencil_op(acc(i + offsets[K])...);
  }

  // 2D
  template<typename A, typename T, std::size_t... K>
  FLECSI_INLINE_TARGET static auto gather(const A & acc,
    const int & axis,
    const T & i,
    const T & j,
    std::index_sequence<K...>) {
    return Derived::stencil_op(shift(acc, axis, offsets[K], i, j)...);
  }

  // 3D
  template<typename A, typename T, std::size_t... K>
  FLECSI_INLINE_TARGET static auto gather(const A & acc,
    const int & axis,
    const T & i,
    const T & j,
    const T & k,
    std::index_sequence<K...>) {
    return Derived::stencil_op(shift(acc, axis, offsets[K], i, j, k)...);
  }

  // 2D
  template<typename A, typename T>
  FLECSI_INLINE_TARGET static auto
  shift(const A & acc, int axis, int o, const T & i, const T & j) {
    return axis == 0 ? acc(i + o, j) : acc(i, j + o);
  }

  // 3D
  template<typename A, typename T>
  FLECSI_INLINE_TARGET static auto
  shift(const A & acc, int axis, int o, const T & i, const T & j, const T & k) {
    return axis == 0 ? acc(i + o, j, k)
                     : (axis == 1 ? acc(i, j + o, k) : acc(i, j, k + o));
  }
};

} // namespace detail

struct genminmod : detail::stencil<genminmod, 3> {

  template<typename T>
  FLECSI_INLINE_TARGET static auto stencil_op(T u_minus, T u_0, T u_plus) {
    using std::abs;
    using std::copysign;
    using std::min;

    const double theta{1.5};
    const T s_0 = 0.5 * (u_plus - u_minus);
    const T s_m = theta * (u_0 - u_minus);
    const T s_p = theta * (u_plus - u_0);

    // Where s_m and s_p share a sign, s_0 shares it too (s_m + s_p =
    // 2*theta*s_0) so 0.5 * (sign(s_m) + sign(s_p)) reproduces the original
    // sign(s_0) factor, and collapses to zero at an extremum. Bit-identical to
    // the branched form.
    const T slope = 0.5 * (copysign(1.0, s_m) + copysign(1.0, s_p)) *
                    min(abs(s_0), min(abs(s_m), abs(s_p)));

    return std::make_tuple(u_0 - 0.5 * slope, u_0 + 0.5 * slope);
  }
}; // genminmod

struct minmod : detail::stencil<minmod, 3> {

  template<typename T>
  FLECSI_INLINE_TARGET static auto stencil_op(T u_minus, T u_0, T u_plus) {
    using std::abs;
    using std::copysign;
    using std::min;

    const T s_m = u_0 - u_minus;
    const T s_p = u_plus - u_0;

    const T slope =
      0.5 * (copysign(1.0, s_p) + copysign(1.0, s_m)) * min(abs(s_p), abs(s_m));

    return std::make_tuple(u_0 - 0.5 * slope, u_0 + 0.5 * slope);
  }
}; // minmod

struct ppm4 : detail::stencil<ppm4, 5> {

  template<typename T>
  FLECSI_INLINE_TARGET static auto
  stencil_op(T u_minus_2, T u_minus, T u_0, T u_plus, T u_plus_2) {

    using std::abs;
    using std::copysign;
    using std::max;
    using std::min;

    const T one_over_twelve(0.08333333333333333);

    // Unlimited 4-th order reconstruction.
    T u_plus_12 =
      (7.0 * (u_0 + u_plus) - (u_plus_2 + u_minus)) * one_over_twelve;
    T u_minus_12 =
      (7.0 * (u_0 + u_minus) - (u_minus_2 + u_plus)) * one_over_twelve;

    // Limit the unlimited reconstruction values to lie between adjacent cell
    // averages.
    u_plus_12 = min(max(u_plus_12, min(u_0, u_plus)), max(u_0, u_plus));
    u_minus_12 = min(max(u_minus_12, min(u_0, u_minus)), max(u_0, u_minus));

    const T delta_u_plus = u_plus_12 - u_0;
    const T delta_u_minus = u_0 - u_minus_12;

    // Monotonization: cap each half-difference at twice the opposite one,
    // keeping its sign (min on magnitudes reproduces the per-branch clamp).
    const T up =
      u_0 +
      copysign(min(abs(delta_u_plus), 2.0 * abs(delta_u_minus)), delta_u_plus);
    const T um =
      u_0 -
      copysign(min(abs(delta_u_minus), 2.0 * abs(delta_u_plus)), delta_u_minus);

    // At a local extremum (opposite-signed differences) flatten to u_0. The
    // mask is 1 where delta_u_plus * delta_u_minus >= 0, else 0.
    const T keep = 0.5 * (1.0 + copysign(1.0, delta_u_plus * delta_u_minus));

    return std::make_tuple(u_0 + keep * (um - u_0), u_0 + keep * (up - u_0));
  }
}; // ppm4

struct weno5z : detail::stencil<weno5z, 5> {

  template<typename T>
  FLECSI_INLINE_TARGET static auto
  stencil_op(T u_minus_2, T u_minus, T u_0, T u_plus, T u_plus_2) {

    const double epsilon = 1.0e-15;

    using std::abs;
    using utils::sqr;

    const std::array<T, 3> beta{
      1.0833333333333333 * sqr(u_minus_2 - 2.0 * u_minus + u_0) +
        0.25 * sqr(u_minus_2 - 4.0 * u_minus + 3.0 * u_0),
      1.0833333333333333 * sqr(u_minus - 2.0 * u_0 + u_plus) +
        0.25 * sqr(u_plus - u_minus),
      1.0833333333333333 * sqr(u_plus_2 - 2.0 * u_plus + u_0) +
        0.25 * sqr(u_plus_2 - 4.0 * u_plus + 3.0 * u_0)};

    const T tau5{abs(beta[2] - beta[0])};

    const std::array<T, 3> epsilon_k{
      epsilon * (1.0 + abs(u_0) + abs(u_minus) + abs(u_minus_2)),
      epsilon * (1.0 + abs(u_0) + abs(u_minus) + abs(u_plus)),
      epsilon * (1.0 + abs(u_0) + abs(u_plus) + abs(u_plus_2))};

    const std::array<T, 3> nw_buffer{1.0 + sqr(tau5 / (beta[0] + epsilon_k[0])),
      1.0 + sqr(tau5 / (beta[1] + epsilon_k[1])),
      1.0 + sqr(tau5 / (beta[2] + epsilon_k[2]))};

    const std::array<T, 3> alpha_upper{
      nw_buffer[0], 10.0 * nw_buffer[1], 5.0 * nw_buffer[2]};
    const std::array<T, 3> alpha_lower{
      nw_buffer[2], 10.0 * nw_buffer[1], 5.0 * nw_buffer[0]};
    const T alpha_norm_upper = alpha_upper[0] + alpha_upper[1] + alpha_upper[2];
    const T alpha_norm_lower = alpha_lower[0] + alpha_lower[1] + alpha_lower[2];

    // reconstruction stencils
    const std::array<T, 3> recons_stencils_upper{
      0.375 * u_minus_2 - 1.25 * u_minus + 1.875 * u_0,
      -0.125 * u_minus + 0.75 * u_0 + 0.375 * u_plus,
      0.375 * u_0 + 0.75 * u_plus - 0.125 * u_plus_2};
    const std::array<T, 3> recons_stencils_lower{
      0.375 * u_plus_2 - 1.25 * u_plus + 1.875 * u_0,
      -0.125 * u_plus + 0.75 * u_0 + 0.375 * u_minus,
      0.375 * u_0 + 0.75 * u_minus - 0.125 * u_minus_2};

    // reconstructed solutions
    return std::make_tuple((alpha_lower[0] * recons_stencils_lower[0] +
                             alpha_lower[1] * recons_stencils_lower[1] +
                             alpha_lower[2] * recons_stencils_lower[2]) /
                             alpha_norm_lower,
      (alpha_upper[0] * recons_stencils_upper[0] +
        alpha_upper[1] * recons_stencils_upper[1] +
        alpha_upper[2] * recons_stencils_upper[2]) /
        alpha_norm_upper);
  }
}; // weno5z

} // namespace spec::reconstructors

#endif // SPEC_RECONSTRUCT_HH
