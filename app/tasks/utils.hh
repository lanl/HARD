#ifndef HARD_TASKS_UTIL_HH
#define HARD_TASKS_UTIL_HH

#include "../state.hh"
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>

#include <sstream>

namespace hard::tasks::util {

/*
  Find the specific internal energy consitent with the given density and
  pressure.
 */
template<typename E>
auto
find_sie(E const & eos,
  double_t r,
  double_t p,
  double_t g = std::numeric_limits<double>::min()) {
  using namespace RootFinding1D;
  auto kernel = [&eos, r](double_t e) { return eos.pRhoSie(r, e); };
  double_t sie{std::numeric_limits<double>::min()};
  const double_t min{eos.eRhoT(r, 1.0e-50)};
  const double_t max{eos.eRhoT(r, 1.0e20)};
  g = g == std::numeric_limits<double>::min() ? 0.5 * (max - min) : g;
  auto s = regula_falsi(kernel, p, g, min, max, 1.0e-12, 1.0e-12, sie);
  flog_assert(s == Status::SUCCESS,
    "speciic internal energy root finder failed for r = " << r
                                                          << " and p = " << p);
  return sie;
} // find_sie

template<dm::domain DM, std::size_t D>
inline void
print_conserved(typename mesh<D>::template accessor<ro> m,
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
    auto rE = m.template mdcolex<is::cells>(rE_a);
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
print_primitives(typename mesh<D>::template accessor<ro> m,
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
} // print_conserved

template<class M, typename IT>
FLECSI_INLINE_TARGET auto
get_mdiota_policy(const M & m, const IT & it1, const IT & it2) {
  using flecsi::exec::mdiota_view;
  using flecsi::exec::sub_range;

  auto b1 = *it1.begin();
  auto e1 = *it1.end();

  auto b2 = *it2.begin();
  auto e2 = *it2.end();

  return mdiota_view(m, sub_range{b1, e1}, sub_range{b2, e2});
} // get_mdiota_policy

template<class M, typename IT>
FLECSI_INLINE_TARGET auto
get_mdiota_policy(const M & m, const IT & it1, const IT & it2, const IT & it3) {
  using flecsi::exec::mdiota_view;
  using flecsi::exec::sub_range;

  auto b1 = *it1.begin();
  auto e1 = *it1.end();

  auto b2 = *it2.begin();
  auto e2 = *it2.end();

  auto b3 = *it3.begin();
  auto e3 = *it3.end();

  return mdiota_view(
    m, sub_range{b1, e1}, sub_range{b2, e2}, sub_range{b3, e3});
} // get_mdiota_policy

} // namespace hard::tasks::util

#endif // HARD_TASKS_UTIL_HH
