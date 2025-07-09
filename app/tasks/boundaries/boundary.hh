#ifndef TASKS_BOUNDARIES_HH
#define TASKS_BOUNDARIES_HH

#include <any>
#include <cstddef>

#include "dirichlet.hh"
#include "flow.hh"
#include "reflective.hh"

namespace hard::tasks {

using hard::tasks::util::bl;

/*----------------------------------------------------------------------------*
  Boundary input translation.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
typename mesh<D>::periodic_axes
init_boundaries(flecsi::exec::cpu,
  typename single<typename mesh<D>::bmap>::template accessor<wo> bmap_a,
  std::array<std::array<bd::boundary_type, 2>, D> bnds) {
  auto & bmap = *bmap_a;
  constexpr auto periodic = mesh<D>::boundary_type::periodic;

  typename mesh<D>::periodic_axes p;
  p[ax::x] = false;
  if constexpr(D == 2 || D == 3) {
    p[ax::y] = false;
  }
  else if constexpr(D == 3) {
    p[ax::z] = false;
  } // if

  bmap[ax::x][bl::low] = bnds[ax::x][bl::low];
  bmap[ax::x][bl::high] = bnds[ax::x][bl::high];

  if((bnds[ax::x][bl::low] == periodic && bnds[ax::x][bl::high] != periodic) ||
     (bnds[ax::x][bl::low] != periodic && bnds[ax::x][bl::high] == periodic))
    flog_fatal("Wrong periodicity in the boundaries.");
  p[ax::x] =
    bnds[ax::x][bl::low] == periodic && bnds[ax::x][bl::high] == periodic;

  if constexpr(D == 2 || D == 3) {
    bmap[ax::y][bl::low] = bnds[ax::y][bl::low];
    bmap[ax::y][bl::high] = bnds[ax::y][bl::high];
    if((bnds[ax::y][bl::low] == periodic &&
         bnds[ax::y][bl::high] != periodic) ||
       (bnds[ax::y][bl::low] != periodic && bnds[ax::y][bl::high] == periodic))
      flog_fatal("Wrong periodicity in the boundaries.");

    p[ax::y] =
      bnds[ax::y][bl::low] == periodic && bnds[ax::y][bl::high] == periodic;
  } // if

  if constexpr(D == 3) {
    bmap[ax::z][bl::low] = bnds[ax::z][bl::low];
    bmap[ax::z][bl::high] = bnds[ax::z][bl::high];
    if((bnds[ax::z][bl::low] == periodic &&
         bnds[ax::z][bl::high] != periodic) ||
       (bnds[ax::z][bl::low] != periodic && bnds[ax::z][bl::high] == periodic))
      p[ax::z] =
        bnds[ax::z][bl::low] == periodic && bnds[ax::z][bl::high] == periodic;
  } // if

  return p;
} // boundaries

template<std::size_t D, typename T>
void
apply_boundary(flecsi::exec::accelerator & s,
  typename mesh<D>::template accessor<ro> & m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> & bmap_a,
  std::vector<typename field<T>::template accessor<rw, ro>> & f_a,
  double value = 0) noexcept {

  using hard::tasks::util::get_mdiota_policy;

  const size_t ghost_zone_size = m.ghost_zone_size();

  const flow<D> f(ghost_zone_size);
  const reflective<D> r(ghost_zone_size);
  const dirichlet<D> d(ghost_zone_size);

  for(auto f_c : f_a) {
    auto f_acc = m.template mdcolex<is::cells>(f_c);
    if constexpr(D == 1) {
      std::array<int, 2> levels = {
        m.template is_low<ax::x>() ? bl::low : bl::none,
        m.template is_high<ax::x>() ? bl::high : bl::none};
      for(auto l : levels)
        if(l != bl::none) {
          std::size_t i = l == bl::low ? 0 : m.template size<ax::x, dm::all>();
          flecsi::util::iota_view policy{0, 1}; // default execution space
          s.executor().forall(x, policy) {
            const typename mesh<D>::bmap & bm = *bmap_a;
            (void)x; // remove compiler warning for unused
            if(bm[0][l] == bd::flow)
              f(f_acc, i, l);
            if(bm[0][l] == bd::reflecting)
              r(f_acc, i, l);
            if(bm[0][l] == bd::dirichlet)
              d(f_acc, value, i, l);
          };
        }
    }
    else if constexpr(D == 2) {
      std::array<int, 2> levels_x = {
        m.template is_low<ax::x>() ? bl::low : bl::none,
        m.template is_high<ax::x>() ? bl::high : bl::none};
      std::array<int, 2> levels_y = {
        m.template is_low<ax::y>() ? bl::low : bl::none,
        m.template is_high<ax::y>() ? bl::high : bl::none};
      for(auto l : levels_x)
        if(l != bl::none) {
          std::size_t i = l == bl::low ? 0 : m.template size<ax::x, dm::all>();
          s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
            const typename mesh<D>::bmap & bm = *bmap_a;
            if(bm[0][l] == bd::flow)
              f(ax::x, f_acc, i, j, l);
            if(bm[0][l] == bd::reflecting)
              r(ax::x, f_acc, i, j, l);
            if(bm[0][l] == bd::dirichlet)
              d(ax::x, f_acc, value, i, j, l);
          };
        }
      for(auto l : levels_y)
        if(l != bl::none) {
          std::size_t j = l == bl::low ? 0 : m.template size<ax::y, dm::all>();
          s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
            const typename mesh<D>::bmap & bm = *bmap_a;
            if(bm[1][l] == bd::flow)
              f(ax::y, f_acc, i, j, l);
            if(bm[1][l] == bd::reflecting)
              r(ax::y, f_acc, i, j, l);
            if(bm[1][l] == bd::dirichlet)
              d(ax::y, f_acc, value, i, j, l);
          }; // forall
        }
    }
    else if constexpr(D == 3) {
      std::array<int, 2> levels_x = {
        m.template is_low<ax::x>() ? bl::low : bl::none,
        m.template is_high<ax::x>() ? bl::high : bl::none};
      std::array<int, 2> levels_y = {
        m.template is_low<ax::y>() ? bl::low : bl::none,
        m.template is_high<ax::y>() ? bl::high : bl::none};
      std::array<int, 2> levels_z = {
        m.template is_low<ax::z>() ? bl::low : bl::none,
        m.template is_high<ax::z>() ? bl::high : bl::none};
      for(auto l : levels_x)
        if(l != bl::none) {
          auto mdpolicy_zy = get_mdiota_policy(f_acc,
            m.template cells<ax::z, dm::quantities>(),
            m.template cells<ax::y, dm::quantities>());
          std::size_t i = l == bl::low ? 0 : m.template size<ax::x, dm::all>();
          s.executor().forall(kj, mdpolicy_zy) {
            const typename mesh<D>::bmap & bm = *bmap_a;
            auto [k, j] = kj;
            if(bm[0][l] == bd::flow)
              f(ax::x, f_acc, i, j, k, l);
            if(bm[0][l] == bd::reflecting)
              r(ax::x, f_acc, i, j, k, l);
            if(bm[0][l] == bd::dirichlet)
              d(ax::x, f_acc, value, i, j, k, l);
          };
        }
      for(auto l : levels_y)
        if(l != bl::none) {
          auto mdpolicy_zx = get_mdiota_policy(f_acc,
            m.template cells<ax::z, dm::quantities>(),
            m.template cells<ax::x, dm::quantities>());
          std::size_t j = l == bl::low ? 0 : m.template size<ax::y, dm::all>();
          s.executor().forall(ki, mdpolicy_zx) {
            const typename mesh<D>::bmap & bm = *bmap_a;
            auto [k, i] = ki;
            if(bm[1][l] == bd::flow)
              f(ax::y, f_acc, i, j, k, l);
            if(bm[1][l] == bd::reflecting)
              r(ax::y, f_acc, i, j, k, l);
            if(bm[1][l] == bd::dirichlet)
              d(ax::y, f_acc, value, i, j, k, l);
          };
        }
      for(auto l : levels_z)
        if(l != bl::none) {
          auto mdpolicy_yx = get_mdiota_policy(f_acc,
            m.template cells<ax::y, dm::quantities>(),
            m.template cells<ax::x, dm::quantities>());
          std::size_t k = l == bl::low ? 0 : m.template size<ax::z, dm::all>();
          s.executor().forall(ji, mdpolicy_yx) {
            const typename mesh<D>::bmap & bm = *bmap_a;
            auto [j, i] = ji;
            if(bm[2][l] == bd::flow)
              f(ax::z, f_acc, i, j, k, l);
            if(bm[2][l] == bd::reflecting)
              r(ax::z, f_acc, i, j, k, l);
            if(bm[2][l] == bd::dirichlet)
              d(ax::z, f_acc, value, i, j, k, l);
          };
        }
    }
  }
}

template<std::size_t D>
void
apply_boundaries(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> bmap_a,
  std::vector<field<double>::accessor<rw, ro>> f_a,
  std::vector<typename field<vec<D>>::template accessor<rw, ro>>
    fv_a) noexcept {
  apply_boundary<D, double>(s, m, bmap_a, f_a);
  apply_boundary<D, vec<D>>(s, m, bmap_a, fv_a);
}

template<std::size_t D>
void
apply_boundaries_scalar(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> bmap_a,
  std::vector<field<double>::accessor<rw, ro>> f_a) noexcept {
  apply_boundary<D, double>(s, m, bmap_a, f_a);
}

template<std::size_t D>
void
apply_boundaries_vector(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> bmap_a,
  std::vector<typename field<vec<D>>::template accessor<rw, ro>> f_a) noexcept {
  apply_boundary<D, vec<D>>(s, m, bmap_a, f_a);
}

template<std::size_t D>
void
apply_dirichlet_boundaries(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> bmap_a,
  std::vector<field<double>::accessor<rw, ro>> f_a,
  single<double>::accessor<ro> value) noexcept {
  apply_boundary<D, double>(s, m, bmap_a, f_a, value);
}

} // namespace hard::tasks
#endif
