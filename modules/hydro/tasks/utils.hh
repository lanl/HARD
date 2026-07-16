#ifndef HARD_MODULE_HYDRO_TASKS_UTIL_HH
#define HARD_MODULE_HYDRO_TASKS_UTIL_HH

#include "../constants.hh"
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>

#include <sstream>

namespace hard::tasks::util {

enum bl { low, high, none };

/*
  Find the specific internal energy consitent with the given density and
  pressure.
 */
template<typename E>
FLECSI_INLINE_TARGET auto
find_sie(E const & eos,
  const double_t r,
  const double_t p,
  double_t g = std::numeric_limits<double>::min()) {
  using namespace RootFinding1D;
  auto kernel = [&eos, r](double_t e) { return eos.pRhoSie(r, e); };
  double_t sie{std::numeric_limits<double>::min()};
  const double_t min{eos.eRhoT(r, 1.0e-50)};
  const double_t max{eos.eRhoT(r, 1.0e20)};
  g = g == std::numeric_limits<double>::min() ? sqrt((max * min)) : g;
  [[maybe_unused]] auto s =
    regula_falsi(kernel, p, g, min, max, 1.0e-12, 1.0e-12, sie);
  assert(s == Status::SUCCESS && "specific internal energy root finder failed");
  return sie;
} // find_sie

template<dm::domain DM, std::size_t D>
inline void
print_conserved(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  typename field<vec<D>>::template accessor<ro, ro> ru_a,
  field<double>::accessor<ro, ro> rE_a,
  flecsi::util::id zslice) {
  {
    auto r = m.template mdcolex<is::cells>(r_a);
    std::stringstream ss;
    ss << "DENSITY:" << std::endl;
    for(auto j : m.template cells<ax::y, DM, true>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << r(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto ru = m.template mdcolex<is::cells>(ru_a);
    std::stringstream ss;
    ss << "MOMENTUM:" << std::endl;
    for(auto j : m.template cells<ax::y, DM, true>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << ru(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto r_e = m.template mdcolex<is::cells>(rE_a);
    std::stringstream ss;
    ss << "TOTAL ENERGY:" << std::endl;
    for(auto j : m.template cells<ax::y, DM, true>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << rE(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
} // print_conserved

template<dm::domain DM, std::size_t D>
inline void
print_vec_field(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename field<vec<D>>::template accessor<ro, ro> u_a) {
  std::stringstream ss;

  if constexpr(D == 1) {
    auto u = m.template mdcolex<is::cells>(u_a);
    for(auto i : m.template cells<ax::x, DM>()) {
      ss << u(i) << " ";
    } // for
    ss << std::endl;
  }
  else if constexpr(D == 2) {
    auto u = m.template mdcolex<is::cells>(u_a);
    for(auto j : m.template cells<ax::y, DM>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << u(i, j) << " ";
      } // for
      ss << std::endl;
    }
  }

  flog(info) << ss.str() << std::endl;
} // print_vec_field

template<dm::domain DM, std::size_t D>
inline void
print_scal_field(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, ro> u_a) {
  std::stringstream ss;

  if constexpr(D == 1) {
    auto u = m.template mdcolex<is::cells>(u_a);
    for(auto i : m.template cells<ax::x, DM>()) {
      ss << u(i) << " ";
    } // for
    ss << std::endl;
  }
  else if constexpr(D == 2) {
    auto u = m.template mdcolex<is::cells>(u_a);
    for(auto j : m.template cells<ax::y, DM>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << u(i, j) << " ";
      } // for
      ss << std::endl;
    }
  }

  flog(info) << ss.str() << std::endl;
} // print_scal_field

template<dm::domain DM, std::size_t D>
inline void
print_primitives(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  typename field<vec<D>>::template accessor<ro, ro> u_a,
  field<double>::accessor<ro, ro> p_a,
  flecsi::util::id zslice) {
  {
    auto u = m.template mdcolex<is::cells>(u_a);
    std::stringstream ss;
    ss << "VELOCITY:" << std::endl;
    for(auto j : m.template cells<ax::y, DM>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << u(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
  {
    auto p = m.template mdcolex<is::cells>(p_a);
    std::stringstream ss;
    ss << "PRESSURE:" << std::endl;
    for(auto j : m.template cells<ax::y, DM>()) {
      for(auto i : m.template cells<ax::x, DM>()) {
        ss << p(i, j, zslice) << " ";
      } // for
      ss << std::endl;
    } // for

    flog(info) << ss.str() << std::endl;
  }
} // print_primitives

} // namespace hard::tasks::util

#endif // HARD_MODULE_HYDRO_TASKS_UTIL_HH
