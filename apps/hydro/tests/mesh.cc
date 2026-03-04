#include "flecsi/data.hh"
#include "flecsi/execution.hh"
#include "flecsi/topology.hh"
#include "flecsi/util/unit.hh"

#include "../../../modules/spec/mesh.hh"

using namespace flecsi;
using namespace flecsi::data;
using namespace flecsi::exec;

using namespace spec;

template<size_t D>
struct state {
  typename mesh<D>::ptr m;
  static inline const field<int>::definition<mesh<D>> field;
};

template<size_t D>
int
task(typename mesh<D>::template accessor<ro>,
  field<int>::accessor<wo, wo>) noexcept {
  return 0;
}

template<std::size_t D>
typename mesh<D>::periodic_axes
init_boundaries() noexcept {
  typename mesh<D>::periodic_axes p;
  p[ax::x] = false;
  if constexpr(D == 2 || D == 3) {
    p[ax::y] = false;
  }
  else if constexpr(D == 3) {
    p[ax::z] = false;
  } // if
  return p;
}

template<std::size_t D>
void
init_mesh(scheduler & s, typename mesh<D>::ptr & m) {

  std::array<std::array<bd::boundary_type, 2>, D> bnds;
  bnds[ax::x][bd::low] = bd::periodic;
  bnds[ax::x][bd::high] = bd::periodic;
  if(D == 2 || D == 3) {
    bnds[ax::y][bd::low] = bd::periodic;
    bnds[ax::y][bd::high] = bd::periodic;
  }
  if(D == 3) {
    bnds[ax::z][bd::low] = bd::periodic;
    bnds[ax::z][bd::high] = bd::periodic;
  }

  auto bf = execute<init_boundaries<D>>();

  std::size_t levels = 5;

  std::size_t lowest_level = levels;
  std::size_t min_highest_level = levels;
  std::size_t max_num_levels = min_highest_level - lowest_level + 1;

  typename mesh<D>::grect geom;
  geom[0][0] = 0.;
  geom[0][1] = 1.;
  if(D == 2 || D == 3) {
    geom[1][0] = 0.;
    geom[1][1] = 1.;
  }
  if(D == 3) {
    geom[2][0] = 0.;
    geom[2][1] = 1.;
  }

  for(std::size_t i{0}; i < max_num_levels; i++) {
    typename mesh<D>::gcoord axis_extents(D);
    axis_extents[ax::x] = 1 << (levels - i);
    if(D == 2 || D == 3) {
      axis_extents[ax::y] = 1 << (levels - i);
    }
    if(D == 3) {
      axis_extents[ax::z] = 1 << (levels - i);
    } // if

    s.allocate(m,
      typename mesh<D>::mpi_coloring(
        s, s.runtime().processes(), axis_extents, bf.get()),
      geom);
  }
}

int
test_driver(scheduler & s) {
  UNIT() {
    // 1D
    state<1> d1;
    init_mesh<1>(s, d1.m);
    EXPECT_EQ(s.test<task<1>>(*d1.m, d1.field(*d1.m)), 0);
    // 2D
    state<2> d2;
    init_mesh<2>(s, d2.m);
    EXPECT_EQ(s.test<task<2>>(*d2.m, d2.field(*d2.m)), 0);
  };
} // kernel_driver

util::unit::driver<test_driver> driver;
