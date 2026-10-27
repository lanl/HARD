#ifndef HARD_COMMON_TASKS_BOUNDARIES_HH
#define HARD_COMMON_TASKS_BOUNDARIES_HH

#include "../modules/common/tasks/utils.hh"
#include "spec/boundaries.hh"

#include "dirichlet.hh"
#include "flow.hh"
#include "reflective.hh"

namespace hard::tasks {

struct init_boundaries_task : annotation_task_region {
  inline static const std::string name{"init_boundaries"};
};
struct apply_boundaries_task : annotation_task_region {
  inline static const std::string name{"apply_boundaries"};
};
struct apply_boundaries_scalar_task : annotation_task_region {
  inline static const std::string name{"apply_boundaries_scalar"};
};

using hard::tasks::util::bl;

/*----------------------------------------------------------------------------*
  Dirichlet value of type T from the components stored in a spec.
 *----------------------------------------------------------------------------*/

template<typename T>
T
bc_value(std::array<double, 3> const & v) {
  if constexpr(std::is_arithmetic_v<T>)
    return static_cast<T>(v[0]);
  else if constexpr(T::Dim == 1)
    return T(v[0]);
  else if constexpr(T::Dim == 2)
    return T(v[0], v[1]);
  else
    return T(v[0], v[1], v[2]);
}

/*----------------------------------------------------------------------------*
  Apply each field's own boundary conditions. specs[n] belongs to f_a[n].
 *----------------------------------------------------------------------------*/

template<std::size_t D, typename T>
void
apply_boundary(flecsi::exec::accelerator & s,
  typename mesh<D>::template accessor<ro> & m,
  std::vector<typename field<T>::template accessor<rw, ro>> & f_a,
  std::vector<spec::bc_spec<D>> const & specs) noexcept {

  using common::tasks::utils::get_mdiota_policy;

  flog_assert(f_a.size() == specs.size(), "one boundary spec per field");

  const size_t ghost_zone_size = m.ghost_zone_size();

  const flow<D> f(ghost_zone_size);
  const reflective<D> r(ghost_zone_size);
  const dirichlet<D> d(ghost_zone_size);

  for(std::size_t n = 0; n < f_a.size(); ++n) {
    auto f_acc = m.template mdcolex<is::cells>(f_a[n]);
    const spec::bc_spec<D> sp = specs[n];

    if constexpr(D == 1) {
      std::array<int, 2> levels = {
        m.template is_low<ax::x>() ? bl::low : bl::none,
        m.template is_high<ax::x>() ? bl::high : bl::none};
      for(auto l : levels)
        if(l != bl::none) {
          const spec::bc_side bs = sp[ax::x][l];
          const T v = bc_value<T>(bs.value);
          std::size_t i = l == bl::low ? 0 : m.template size<ax::x, dm::all>();
          flecsi::util::iota_view policy{0, 1}; // default execution space
          s.executor().named("boundary_x").forall(x, policy) {
            (void)x; // remove compiler warning for unused
            if(bs.type == bd::flow)
              f(f_acc, i, l);
            if(bs.type == bd::reflecting)
              r(f_acc, i, l);
            if(bs.type == bd::dirichlet)
              d(f_acc, v, i, l);
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
          const spec::bc_side bs = sp[ax::x][l];
          const T v = bc_value<T>(bs.value);
          std::size_t i = l == bl::low ? 0 : m.template size<ax::x, dm::all>();
          s.executor()
            .named("boundary_x")
            .forall(j, (m.template cells<ax::y, dm::quantities>())) {
            if(bs.type == bd::flow)
              f(ax::x, f_acc, i, j, l);
            if(bs.type == bd::reflecting)
              r(ax::x, f_acc, i, j, l);
            if(bs.type == bd::dirichlet)
              d(ax::x, f_acc, v, i, j, l);
          };
        }
      for(auto l : levels_y)
        if(l != bl::none) {
          const spec::bc_side bs = sp[ax::y][l];
          const T v = bc_value<T>(bs.value);
          std::size_t j = l == bl::low ? 0 : m.template size<ax::y, dm::all>();
          s.executor()
            .named("boundary_y")
            .forall(i, (m.template cells<ax::x, dm::quantities>())) {
            if(bs.type == bd::flow)
              f(ax::y, f_acc, i, j, l);
            if(bs.type == bd::reflecting)
              r(ax::y, f_acc, i, j, l);
            if(bs.type == bd::dirichlet)
              d(ax::y, f_acc, v, i, j, l);
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
          const spec::bc_side bs = sp[ax::x][l];
          const T v = bc_value<T>(bs.value);
          auto mdpolicy_zy = get_mdiota_policy(f_acc,
            m.template cells<ax::z, dm::quantities>(),
            m.template cells<ax::y, dm::quantities>());
          std::size_t i = l == bl::low ? 0 : m.template size<ax::x, dm::all>();
          s.executor().named("boundary_x").forall(kj, mdpolicy_zy) {
            auto [k, j] = kj;
            if(bs.type == bd::flow)
              f(ax::x, f_acc, i, j, k, l);
            if(bs.type == bd::reflecting)
              r(ax::x, f_acc, i, j, k, l);
            if(bs.type == bd::dirichlet)
              d(ax::x, f_acc, v, i, j, k, l);
          };
        }
      for(auto l : levels_y)
        if(l != bl::none) {
          const spec::bc_side bs = sp[ax::y][l];
          const T v = bc_value<T>(bs.value);
          auto mdpolicy_zx = get_mdiota_policy(f_acc,
            m.template cells<ax::z, dm::quantities>(),
            m.template cells<ax::x, dm::quantities>());
          std::size_t j = l == bl::low ? 0 : m.template size<ax::y, dm::all>();
          s.executor().named("boundary_y").forall(ki, mdpolicy_zx) {
            auto [k, i] = ki;
            if(bs.type == bd::flow)
              f(ax::y, f_acc, i, j, k, l);
            if(bs.type == bd::reflecting)
              r(ax::y, f_acc, i, j, k, l);
            if(bs.type == bd::dirichlet)
              d(ax::y, f_acc, v, i, j, k, l);
          };
        }
      for(auto l : levels_z)
        if(l != bl::none) {
          const spec::bc_side bs = sp[ax::z][l];
          const T v = bc_value<T>(bs.value);
          auto mdpolicy_yx = get_mdiota_policy(f_acc,
            m.template cells<ax::y, dm::quantities>(),
            m.template cells<ax::x, dm::quantities>());
          std::size_t k = l == bl::low ? 0 : m.template size<ax::z, dm::all>();
          s.executor().named("boundary_z").forall(ji, mdpolicy_yx) {
            auto [j, i] = ji;
            if(bs.type == bd::flow)
              f(ax::z, f_acc, i, j, k, l);
            if(bs.type == bd::reflecting)
              r(ax::z, f_acc, i, j, k, l);
            if(bs.type == bd::dirichlet)
              d(ax::z, f_acc, v, i, j, k, l);
          };
        }
    }
  }
}

/*----------------------------------------------------------------------------*
  Task: scalar fields with their specs, vector fields with their specs.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
void
apply_boundaries(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  std::vector<field<double>::accessor<rw, ro>> f_a,
  std::vector<spec::bc_spec<D>> f_specs,
  std::vector<typename field<vec<D>>::template accessor<rw, ro>> fv_a,
  std::vector<spec::bc_spec<D>> fv_specs,
  std::vector<field<int>::accessor<rw, ro>> fi_a,
  std::vector<spec::bc_spec<D>> fi_specs) noexcept {
  apply_boundary<D, double>(s, m, f_a, f_specs);
  apply_boundary<D, vec<D>>(s, m, fv_a, fv_specs);
  apply_boundary<D, int>(s, m, fi_a, fi_specs);
}

} // namespace hard::tasks

/*----------------------------------------------------------------------------*
  Host side: collect fields and their specs, launch one task.
 *----------------------------------------------------------------------------*/

namespace hard::actions {

template<std::size_t D>
struct boundary_batch {
  using scalar_ref = field<double>::Reference<mesh<D>, is::cells>;
  using vector_ref =
    typename field<vec<D>>::template Reference<mesh<D>, is::cells>;
  using int_ref = field<int>::Reference<mesh<D>, is::cells>;

  explicit boundary_batch(spec::boundaries<D> const & bcs) : bcs_(bcs) {}

  // Configurable field: its spec from the input file, by name.
  template<typename R>
  void add(std::string const & name, R r) {
    push(r, bcs_.spec(name));
  }
  template<typename F, typename M>
  void add(F const & fld, M & m) {
    add(fld.name, fld(m));
  }
  template<typename Tuple, typename M>
  void add_tagged(Tuple const & fields, tag::field_tags tg, M & m) {
    for_each_field(fields, tg, [&](auto const & fld) { add(fld, m); });
  }

  // Non-configurable field: flow on every side (level sets, material ids).
  template<typename R>
  void add_flow(R r) {
    push(r, spec::boundaries<D>::uniform(bd::flow));
  }

  template<typename M>
  void execute(flecsi::scheduler & sc, M & m) const {
    annotation_rguard<tasks::apply_boundaries_task> guard;
    sc.execute<tasks::apply_boundaries<D>>(
      flecsi::exec::on, m, sa_, ss_, va_, vs_, ia_, is_);
  }

private:
  void push(scalar_ref r, spec::bc_spec<D> const & sp) {
    sa_.push_back(r);
    ss_.push_back(sp);
  }
  void push(vector_ref r, spec::bc_spec<D> const & sp) {
    va_.push_back(r);
    vs_.push_back(sp);
  }
  void push(int_ref r, spec::bc_spec<D> const & sp) {
    ia_.push_back(r);
    is_.push_back(sp);
  }

  spec::boundaries<D> const & bcs_;
  std::vector<scalar_ref> sa_;
  std::vector<spec::bc_spec<D>> ss_;
  std::vector<vector_ref> va_;
  std::vector<spec::bc_spec<D>> vs_;
  std::vector<int_ref> ia_;
  std::vector<spec::bc_spec<D>> is_;
};

// Apply the configured boundaries to every field of `s.fields()` carrying
// the tag `tg`, plus flow on each extra field reference in `flow_refs`
// (level sets, material ids), in a single task launch.
template<std::size_t D, typename S, typename... R>
void
apply_tagged_boundaries(flecsi::scheduler & sc,
  S & s,
  tag::field_tags tg,
  R... flow_refs) {
  boundary_batch<D> bb(s.bcs);
  bb.add_tagged(s.fields(), tg, *s.m);
  (bb.add_flow(flow_refs), ...);
  bb.execute(sc, *s.m);
}

// Flow on every side for each field reference, regardless of the input
// file, in a single task launch.
template<std::size_t D, typename S, typename... R>
void
apply_flow_boundaries(flecsi::scheduler & sc, S & s, R... flow_refs) {
  boundary_batch<D> bb(s.bcs);
  (bb.add_flow(flow_refs), ...);
  bb.execute(sc, *s.m);
}

} // namespace hard::actions

#endif // HARD_COMMON_TASKS_BOUNDARIES_HH
