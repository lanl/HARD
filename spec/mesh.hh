#ifndef SPEC_MESH_HH
#define SPEC_MESH_HH

#include "labels.hh"

#include <flecsi/data.hh>
#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include <ranges>

namespace spec {

namespace detail {

template<std::size_t D>
struct help {};

template<>
struct help<1> {
  using axes = flecsi::topo::help::has<ax::x>;
  struct meta_data {
    double xdelta;
  };
};

template<>
struct help<2> {
  using axes = flecsi::topo::help::has<ax::x, ax::y>;
  struct meta_data {
    double xdelta;
    double ydelta;
  };
};

template<>
struct help<3> {
  using axes = flecsi::topo::help::has<ax::x, ax::y, ax::z>;
  struct meta_data {
    double xdelta;
    double ydelta;
    double zdelta;
  };
};
} // namespace detail

/*!
  FIXME
 */
template<std::size_t D>
struct mesh : flecsi::topo::specialization<flecsi::topo::narray, mesh<D>> {

  static constexpr flecsi::util::id ghost_zone_size_ = 3;

  static_assert(D == 1 || D == 2 || D == 3, "unsupported dimension");

  /*--------------------------------------------------------------------------*
    Policy Information.
   *--------------------------------------------------------------------------*/

  static constexpr flecsi::Dimension dimension = D;
  using index_space = is::index_space;
  using index_spaces = flecsi::topo::help::has<is::cells>;
  using axis = ax::axis;
  using axes = typename detail::help<D>::axes;
  using boundary = bd::boundary;
  using boundary_type = bd::boundary_type;
  using coord = typename mesh::base::coord;
  using gcoord = typename mesh::base::gcoord;
  using periodic_axes = std::array<bool, D>;
  using bmap = std::array<std::array<bd::boundary_type, 2 /* low, high*/>, D>;
  using grect = std::array<std::array<double, 2 /* low, high */>, D>;
  using colors = typename mesh::base::colors;
  using coloring = typename mesh::base::coloring;
  using color_coord = std::array<flecsi::Color, D>;
  using axis_definition = typename mesh::base::axis_definition;
  using index_definition = typename mesh::base::index_definition;

  using meta_data = typename detail::help<D>::meta_data;

  template<auto>
  static constexpr std::size_t privilege_count = 2;

  template<auto S, typename T>
  FLECSI_INLINE_TARGET static auto make_ids(T && t) {
    return flecsi::util::transform_view(std::forward<T>(t),
      [](auto const & i) { return flecsi::topo::id<S>(i); });
  } // make_ids

  /*--------------------------------------------------------------------------*
    Interface.
   *--------------------------------------------------------------------------*/

  /// Mesh Interface.
  template<class B>
  struct interface : B {

    FLECSI_INLINE_TARGET std::size_t ghost_zone_size() const {
      return ghost_zone_size_;
    } // global_id

    /// Return axis info for the given axis.
    /// @tparam A The mesh axis.
    template<mesh::axis A>
    FLECSI_INLINE_TARGET typename mesh::base::axis_info axis() const {
      return B::template axis<is::cells, A>();
    }

    /// Return the size for the given axis and domain.
    /// @tparam A  The mesh axis.
    /// @tparam DM The mesh domain.
    template<ax::axis A, dm::domain DM = dm::quantities>
    auto size() const {
      const typename mesh::base::axis_info a = axis<A>();
      if constexpr(DM == dm::quantities) {
        return a.logical;
      }
      else if constexpr(DM == dm::predictor) {
        return a.logical + 2;
      }
      else if constexpr(DM == dm::corrector) {
        return a.logical + 1;
      }
      if constexpr(DM == dm::interior) {
        const bool low = axis<A>().low();
        const bool high = axis<A>().high();

        if(low && high) {
          return a.logical - 2;
        }
        else if(low || high) {
          return a.logical - 1;
        }
        else { /* interior */
          return a.logical;
        }
      }
      else if constexpr(DM == dm::all) {
        return a.layout.extent();
      }
      else if constexpr(DM == dm::global) {
        return a.axis.extent;
      } // if
    } // size

    /// Return a range over the given axis and domain.
    /// @tparam A  The mesh axis.
    /// @tparam DM The mesh domain (excluding \em global).
    ///
    /// The range can be used to iterate over the cells in the given domain,
    /// e.g.:
    /// @code
    ///   for(auto k: m.cells<ax::z, mesh::quantities>()) {
    ///     for(auto j: m.cells<ax::y, mesh::quantities>()) {
    ///       for(auto i: m.cells<ax::x, mesh::quantities>()) {
    ///         u[k][j][i] = 1.0;
    ///       } // for
    ///     } // for
    ///   } // for
    /// @endcode
    template<ax::axis A, dm::domain DM = dm::quantities>
    FLECSI_INLINE_TARGET auto cells() const {
      flecsi::util::id b, e = axis<A>().layout.extent();

      if constexpr(DM == dm::quantities) {
        b = ghost_zone_size_;
        e -= ghost_zone_size_;
      }
      else if constexpr(DM == dm::predictor) {
        b = ghost_zone_size_ - 1;
        e -= (ghost_zone_size_ - 1);
      }
      else if constexpr(DM == dm::corrector) {
        b = ghost_zone_size_;
        e -= (ghost_zone_size_ - 1);
      }
      else if constexpr(DM == dm::interior) {
        b = ghost_zone_size_ + 1;
        e -= ghost_zone_size_ + 1;
      }
      else if constexpr(DM == dm::all) {
        b = 0;
      }
      else if(DM == dm::global) {
        flog_fatal(
          "illegal domain: you cannot iterate over the global domain.");
      } // if

      return make_ids<is::cells>(flecsi::util::iota_view{b, e});

    } // cells

    template<typename... CC>
    void set_geometry(CC... cc) {
      this->policy_meta() = {
        cc...,
      };
    } // set_geometry

    /// Return the global id of \em i for the given axis.
    /// @tparam A The coordinate axis.
    template<ax::axis A>
    FLECSI_INLINE_TARGET flecsi::util::gid global_id(std::size_t i) const {
      return axis<A>().global_id(i);
    } // global_id

    /// Return the mesh spacing for the given axis.
    /// @tparam A  The coordinate axis.
    template<ax::axis A>
    FLECSI_INLINE_TARGET double delta() const {
      if constexpr(A == ax::x) {
        return this->policy_meta().xdelta;
      }
      else if constexpr(D == 2 || (D == 3 && A == ax::y)) {
        return this->policy_meta().ydelta;
      }
      else if constexpr(D == 3 && A == ax::z) {
        return this->policy_meta().zdelta;
      } // if
    } // delta

    /// Return the cell head for the given axis and id. The head is the trailing
    /// interface of the cell.
    /// @tparam A  The coordinate axis.
    template<ax::axis A>
    FLECSI_INLINE_TARGET double head(std::size_t i) const {
      return center<A>(i) - 0.5 * delta<A>();
    } // center

    /// Return the cell center for the given axis and id.
    /// @tparam A  The coordinate axis.
    template<ax::axis A>
    FLECSI_INLINE_TARGET double center(std::size_t i) const {
      return delta<A>() * global_id<A>(i) + 0.5 * delta<A>();
    } // center

    /// Return the cell tail for the given axis and id. The tail is the leading
    /// interface of the cell.
    /// @tparam A  The coordinate axis.
    template<ax::axis A>
    FLECSI_INLINE_TARGET double tail(std::size_t i) const {
      return center<A>(i) + 0.5 * delta<A>();
    } // center

    /// Return true if the current color is a low-edge partition for the given
    /// axis.
    /// @tparam A  The coordinate axis.
    template<ax::axis A>
    bool is_low() {
      return axis<A>().low();
    } // is_low

    /// Return true if the current color is a high-edge partition for the given
    /// axis.
    /// @tparam A  The coordinate axis.
    template<ax::axis A>
    bool is_high() {
      return axis<A>().high();
    } // is_high

    color_coord color_indeces() const {
      if constexpr(D == 1) {
        return {axis<ax::x>().color};
      }
      else if constexpr(D == 2) {
        return {axis<ax::x>().color, axis<ax::y>().color};
      }
      else /* D == 3 */ {
        return {axis<ax::x>().color, axis<ax::y>().color, axis<ax::z>().color};
      } // if
    } // color_indeces

    color_coord axis_colors() const {
      if constexpr(D == 1) {
        return {axis<ax::x>().axis.colors};
      }
      else if constexpr(D == 2) {
        return {axis<ax::x>().axis.colors, axis<ax::y>().axis.colors};
      }
      else /* D == 3 */ {
        return {axis<ax::x>().axis.colors,
          axis<ax::y>().axis.colors,
          axis<ax::z>().axis.colors};
      } // if
    } // axis_colors

  }; // interface

  /*--------------------------------------------------------------------------*
    Color Task.
   *--------------------------------------------------------------------------*/

  template<typename T>
  index_definition
  index_colors(T & num, gcoord & axis_extents, periodic_axes & p) {
    index_definition idef;
    idef.axes = mesh::base::make_axes(num, axis_extents);
    std::size_t ai{0};

    for(auto & a : idef.axes) {
      a.hdepth = ghost_zone_size_;
      a.bdepth = ghost_zone_size_;
      a.periodic = p[ai++];
    } // for

    return {{idef}};
  }

  static coloring
  color(flecsi::Color num_colors, gcoord axis_extents, periodic_axes p) {
    return {{index_colors(num_colors, axis_extents, p)}};
  }

  static coloring color(mesh::base::colors const & color_distribution,
    gcoord axis_extents,
    periodic_axes p) {
    return {{index_colors(color_distribution, axis_extents, p)}};
  } // color

  /*--------------------------------------------------------------------------*
    Initialization.
   *--------------------------------------------------------------------------*/

  static void set_geometry(mesh::template accessor<flecsi::rw> m,
    grect const & g) {
    if constexpr(D == 1) {
      m.set_geometry(
        (g[0][1] - g[0][0]) / m.template size<ax::x, dm::global>());
    }
    else if constexpr(D == 2) {
      m.set_geometry((g[0][1] - g[0][0]) / m.template size<ax::x, dm::global>(),
        (g[1][1] - g[1][0]) / m.template size<ax::y, dm::global>());
    }
    else /* D == 3 */ {
      m.set_geometry((g[0][1] - g[0][0]) / m.template size<ax::x, dm::global>(),
        (g[1][1] - g[1][0]) / m.template size<ax::y, dm::global>(),
        (g[2][1] - g[2][0]) / m.template size<ax::z, dm::global>());
    }
  } // set_geometry

  static void initialize(flecsi::data::topology_slot<mesh> & s,
    coloring const &,
    grect const & geometry) {
    flecsi::execute<set_geometry>(s, geometry);
  } // initialize

}; // struct mesh

} // namespace spec

#endif // SPEC_MESH_HH
