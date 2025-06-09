
#pragma once

#include "../tasks/utils.hh"
#include "../types.hh"
#include <cstddef>

#define LOW 0
#define HIGH 1

namespace hard::tasks {

template<std::size_t D>
void
apply_boundaries(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> bmap_a,
  // primitive variables
  field<double>::accessor<rw, na> r_a,
  typename field<vec<D>>::template accessor<rw, ro> v_a,
  field<double>::accessor<rw, ro> p_a,
  field<double>::accessor<rw, ro> sie_a,
  field<double>::accessor<rw, ro> Erad_a,
  // conservative variables
  typename field<vec<D>>::template accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a) {
  using hard::tasks::util::get_mdiota_policy;

  auto r = m.template mdcolex<is::cells>(r_a);
  auto v = m.template mdcolex<is::cells>(v_a);
  auto p = m.template mdcolex<is::cells>(p_a);
  auto sie = m.template mdcolex<is::cells>(sie_a);
  auto Erad = m.template mdcolex<is::cells>(Erad_a);

  auto ru = m.template mdcolex<is::cells>(ru_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);

  const size_t ghost_zone_size = m.ghost_zone_size();

  if(m.template is_low<ax::x>()) { // NOTE: true only for boundary colors
    if constexpr(D == 1) {
      flecsi::util::iota_view policy{0, 1}; // default execution space
      s.executor().forall(x, policy) {
        (void)x; // remove compiler warning for unused

        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xlow = bm[0][LOW];
        if(xlow == bd::flow) { // NOTE: Same for all threads/warps in a color.
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m) = r(ghost_zone_size);
            v(m) = v(ghost_zone_size);
            p(m) = p(ghost_zone_size);
            sie(m) = sie(ghost_zone_size);

            ru(m) = ru(ghost_zone_size);
            rE(m) = rE(ghost_zone_size);
            Erad(m) = Erad(ghost_zone_size);
          }
        }
        else if(xlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m) = r(2 * ghost_zone_size - 1 - m);
            v(m).x = -v(2 * ghost_zone_size - 1 - m).x;
            p(m) = p(2 * ghost_zone_size - 1 - m);
            sie(m) = sie(2 * ghost_zone_size - 1 - m);

            ru(m).x = -ru(2 * ghost_zone_size - 1 - m).x;
            rE(m) = rE(2 * ghost_zone_size - 1 - m);
            Erad(m) = Erad(2 * ghost_zone_size - 1 - m);
          }
        }
      };
    }
    else if constexpr(D == 2) {
      s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xlow = bm[0][LOW];
        if(xlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j) = r(ghost_zone_size, j);
            v(m, j) = v(ghost_zone_size, j);
            p(m, j) = p(ghost_zone_size, j);
            sie(m, j) = sie(ghost_zone_size, j);

            ru(m, j) = ru(ghost_zone_size, j);
            rE(m, j) = rE(ghost_zone_size, j);
            Erad(m, j) = Erad(ghost_zone_size, j);
          }
        }
        else if(xlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j) = r(2 * ghost_zone_size - 1 - m, j);
            v(m, j).x = -v(2 * ghost_zone_size - 1 - m, j).x;
            v(m, j).y = v(2 * ghost_zone_size - 1 - m, j).y;
            p(m, j) = p(2 * ghost_zone_size - 1 - m, j);
            sie(m, j) = sie(2 * ghost_zone_size - 1 - m, j);

            ru(m, j).x = -ru(2 * ghost_zone_size - 1 - m, j).x;
            ru(m, j).y = ru(2 * ghost_zone_size - 1 - m, j).y;
            rE(m, j) = rE(2 * ghost_zone_size - 1 - m, j);
            Erad(m, j) = Erad(2 * ghost_zone_size - 1 - m, j);
          }
        }
      }; // for
    }
    else /* D == 3 */ {
      auto mdpolicy_zy = get_mdiota_policy(r,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>());

      s.executor().forall(kj, mdpolicy_zy) {
        auto [k, j] = kj;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xlow = bm[0][LOW];
        if(xlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j, k) = r(ghost_zone_size, j, k);
            v(m, j, k) = v(ghost_zone_size, j, k);
            p(m, j, k) = p(ghost_zone_size, j, k);
            sie(m, j, k) = sie(ghost_zone_size, j, k);

            ru(m, j, k) = ru(ghost_zone_size, j, k);
            rE(m, j, k) = rE(ghost_zone_size, j, k);
            Erad(m, j, k) = Erad(ghost_zone_size, j, k);
          }
        }
        else if(xlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j, k) = r(2 * ghost_zone_size - 1 - m, j, k);
            v(m, j, k).x = -v(2 * ghost_zone_size - 1 - m, j, k).x;
            v(m, j, k).y = v(2 * ghost_zone_size - 1 - m, j, k).y;
            v(m, j, k).z = v(2 * ghost_zone_size - 1 - m, j, k).z;
            p(m, j, k) = p(2 * ghost_zone_size - 1 - m, j, k);
            sie(m, j, k) = sie(2 * ghost_zone_size - 1 - m, j, k);

            ru(m, j, k).x = -ru(2 * ghost_zone_size - 1 - m, j, k).x;
            ru(m, j, k).y = ru(2 * ghost_zone_size - 1 - m, j, k).y;
            ru(m, j, k).z = ru(2 * ghost_zone_size - 1 - m, j, k).z;
            rE(m, j, k) = rE(2 * ghost_zone_size - 1 - m, j, k);
            Erad(m, j, k) = Erad(2 * ghost_zone_size - 1 - m, j, k);
          }
        }
      }; // forall
    } // if
  }
  if(m.template is_high<ax::x>()) {
    const std::size_t i = m.template size<ax::x, dm::all>();
    if constexpr(D == 1) {
      flecsi::util::iota_view policy{0, 1}; // default execution space
      s.executor().forall(x, policy) {
        (void)x; // remove compiler warning for unused

        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xhigh = bm[0][HIGH];
        if(xhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m) = r(i - 1 - ghost_zone_size);
            v(i - 1 - m) = v(i - 1 - ghost_zone_size);
            p(i - 1 - m) = p(i - 1 - ghost_zone_size);
            sie(i - 1 - m) = sie(i - 1 - ghost_zone_size);

            ru(i - 1 - m) = ru(i - 1 - ghost_zone_size);
            rE(i - 1 - m) = rE(i - 1 - ghost_zone_size);
            Erad(i - 1 - m) = Erad(i - 1 - ghost_zone_size);
          }
        }
        else if(xhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m) = r(i - 2 * ghost_zone_size + m);
            v(i - 1 - m) = v(i - 2 * ghost_zone_size + m);
            p(i - 1 - m) = p(i - 2 * ghost_zone_size + m);
            sie(i - 1 - m) = sie(i - 2 * ghost_zone_size + m);

            ru(i - 1 - m) = ru(i - 2 * ghost_zone_size + m);
            rE(i - 1 - m) = rE(i - 2 * ghost_zone_size + m);
            Erad(i - 1 - m) = Erad(i - 2 * ghost_zone_size + m);
          }
        }
      };
    }
    else if constexpr(D == 2) {
      s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xhigh = bm[0][HIGH];
        if(xhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j) = r(i - 1 - ghost_zone_size, j);
            v(i - 1 - m, j) = v(i - 1 - ghost_zone_size, j);
            p(i - 1 - m, j) = p(i - 1 - ghost_zone_size, j);
            sie(i - 1 - m, j) = sie(i - 1 - ghost_zone_size, j);

            ru(i - 1 - m, j) = ru(i - 1 - ghost_zone_size, j);
            rE(i - 1 - m, j) = rE(i - 1 - ghost_zone_size, j);
            Erad(i - 1 - m, j) = Erad(i - 1 - ghost_zone_size, j);
          }
        }
        else if(xhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j) = r(i - 2 * ghost_zone_size + m, j);
            v(i - 1 - m, j).x = -v(i - 2 * ghost_zone_size + m, j).x;
            v(i - 1 - m, j).y = v(i - 2 * ghost_zone_size + m, j).y;
            p(i - 1 - m, j) = p(i - 2 * ghost_zone_size + m, j);
            sie(i - 1 - m, j) = sie(i - 2 * ghost_zone_size + m, j);

            ru(i - 1 - m, j).x = -ru(i - 2 * ghost_zone_size + m, j).x;
            ru(i - 1 - m, j).y = ru(i - 2 * ghost_zone_size + m, j).y;
            rE(i - 1 - m, j) = rE(i - 2 * ghost_zone_size + m, j);
            Erad(i - 1 - m, j) = Erad(i - 2 * ghost_zone_size + m, j);
          }
        }
      }; // forall
    }
    else /* D == 3 */ {
      auto mdpolicy_kj = get_mdiota_policy(r,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>());

      s.executor().forall(kj, mdpolicy_kj) {
        auto [k, j] = kj;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xhigh = bm[0][HIGH];
        if(xhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j, k) = r(i - 1 - ghost_zone_size, j, k);
            v(i - 1 - m, j, k) = v(i - 1 - ghost_zone_size, j, k);
            p(i - 1 - m, j, k) = p(i - 1 - ghost_zone_size, j, k);
            sie(i - 1 - m, j, k) = sie(i - 1 - ghost_zone_size, j, k);

            ru(i - 1 - m, j, k) = ru(i - 1 - ghost_zone_size, j, k);
            rE(i - 1 - m, j, k) = rE(i - 1 - ghost_zone_size, j, k);
            Erad(i - 1 - m, j, k) = Erad(i - 1 - ghost_zone_size, j, k);
          }
        }
        else if(xhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j, k) = r(i - 2 * ghost_zone_size + m, j, k);
            v(i - 1 - m, j, k).x = -v(i - 2 * ghost_zone_size + m, j, k).x;
            v(i - 1 - m, j, k).y = v(i - 2 * ghost_zone_size + m, j, k).y;
            v(i - 1 - m, j, k).z = v(i - 2 * ghost_zone_size + m, j, k).z;
            p(i - 1 - m, j, k) = p(i - 2 * ghost_zone_size + m, j, k);
            sie(i - 1 - m, j, k) = sie(i - 2 * ghost_zone_size + m, j, k);

            ru(i - 1 - m, j, k).x = -ru(i - 2 * ghost_zone_size + m, j, k).x;
            ru(i - 1 - m, j, k).y = ru(i - 2 * ghost_zone_size + m, j, k).y;
            ru(i - 1 - m, j, k).z = ru(i - 2 * ghost_zone_size + m, j, k).z;
            rE(i - 1 - m, j, k) = rE(i - 2 * ghost_zone_size + m, j, k);
            Erad(i - 1 - m, j, k) = Erad(i - 2 * ghost_zone_size + m, j, k);
          }
        }

      }; // forall
    } // if
  } // if

  if constexpr(D == 2 || D == 3) {
    if(m.template is_low<ax::y>()) {
      if constexpr(D == 2) {
        s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto ylow = bm[1][LOW];
          if(ylow == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m) = r(i, ghost_zone_size);
              v(i, m) = v(i, ghost_zone_size);
              p(i, m) = p(i, ghost_zone_size);
              sie(i, m) = sie(i, ghost_zone_size);

              ru(i, m) = ru(i, ghost_zone_size);
              rE(i, m) = rE(i, ghost_zone_size);
              Erad(i, m) = Erad(i, ghost_zone_size);
            }
          }
          else if(ylow == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m) = r(i, 2 * ghost_zone_size - 1 - m);
              v(i, m).x = v(i, 2 * ghost_zone_size - 1 - m).x;
              v(i, m).y = -v(i, 2 * ghost_zone_size - 1 - m).y;
              p(i, m) = p(i, 2 * ghost_zone_size - 1 - m);
              sie(i, m) = sie(i, 2 * ghost_zone_size - 1 - m);

              ru(i, m).x = ru(i, 2 * ghost_zone_size - 1 - m).x;
              ru(i, m).y = -ru(i, 2 * ghost_zone_size - 1 - m).y;
              rE(i, m) = rE(i, 2 * ghost_zone_size - 1 - m);
              Erad(i, m) = Erad(i, 2 * ghost_zone_size - 1 - m);
            }
          }
        }; // forall
      }
      else /* D == 3 */ {
        auto mdpolicy_zx = get_mdiota_policy(r,
          m.template cells<ax::z, dm::quantities>(),
          m.template cells<ax::x, dm::quantities>());
        s.executor().forall(ki, mdpolicy_zx) {
          auto [k, i] = ki;

          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto ylow = bm[1][LOW];
          if(ylow == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m, k) = r(i, ghost_zone_size, k);
              v(i, m, k) = v(i, ghost_zone_size, k);
              p(i, m, k) = p(i, ghost_zone_size, k);
              sie(i, m, k) = sie(i, ghost_zone_size, k);

              ru(i, m, k) = ru(i, ghost_zone_size, k);
              rE(i, m, k) = rE(i, ghost_zone_size, k);
              Erad(i, m, k) = Erad(i, ghost_zone_size, k);
            }
          }
          else if(ylow == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m, k) = r(i, 2 * ghost_zone_size - 1 - m, k);
              v(i, m, k).x = v(i, 2 * ghost_zone_size - 1 - m, k).x;
              v(i, m, k).y = -v(i, 2 * ghost_zone_size - 1 - m, k).y;
              v(i, m, k).z = v(i, 2 * ghost_zone_size - 1 - m, k).z;
              p(i, m, k) = p(i, 2 * ghost_zone_size - 1 - m, k);
              sie(i, m, k) = sie(i, 2 * ghost_zone_size - 1 - m, k);

              ru(i, m, k).x = ru(i, 2 * ghost_zone_size - 1 - m, k).x;
              ru(i, m, k).y = -ru(i, 2 * ghost_zone_size - 1 - m, k).y;
              ru(i, m, k).z = ru(i, 2 * ghost_zone_size - 1 - m, k).z;
              rE(i, m, k) = rE(i, 2 * ghost_zone_size - 1 - m, k);
              Erad(i, m, k) = Erad(i, 2 * ghost_zone_size - 1 - m, k);
            }
          }
        }; // forall
      } // if
    }
    if(m.template is_high<ax::y>()) {
      const std::size_t j = m.template size<ax::y, dm::all>();
      if constexpr(D == 2) {
        s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto yhigh = bm[1][HIGH];
          if(yhigh == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m) = r(i, j - 1 - ghost_zone_size);
              v(i, j - 1 - m) = v(i, j - 1 - ghost_zone_size);
              p(i, j - 1 - m) = p(i, j - 1 - ghost_zone_size);
              sie(i, j - 1 - m) = sie(i, j - 1 - ghost_zone_size);

              ru(i, j - 1 - m) = ru(i, j - 1 - ghost_zone_size);
              rE(i, j - 1 - m) = rE(i, j - 1 - ghost_zone_size);
              Erad(i, j - 1 - m) = Erad(i, j - 1 - ghost_zone_size);
            }
          }
          else if(yhigh == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m) = r(i, j - 2 * ghost_zone_size + m);
              v(i, j - 1 - m).x = v(i, j - 2 * ghost_zone_size + m).x;
              v(i, j - 1 - m).y = -v(i, j - 2 * ghost_zone_size + m).y;
              p(i, j - 1 - m) = p(i, j - 2 * ghost_zone_size + m);
              sie(i, j - 1 - m) = sie(i, j - 2 * ghost_zone_size + m);

              ru(i, j - 1 - m).x = ru(i, j - 2 * ghost_zone_size + m).x;
              ru(i, j - 1 - m).y = -ru(i, j - 2 * ghost_zone_size + m).y;
              rE(i, j - 1 - m) = rE(i, j - 2 * ghost_zone_size + m);
              Erad(i, j - 1 - m) = Erad(i, j - 2 * ghost_zone_size + m);
            }
          }
        }; // forall
      }
      else /* D == 3 */ {
        auto mdpolicy_zx = get_mdiota_policy(r,
          m.template cells<ax::z, dm::quantities>(),
          m.template cells<ax::x, dm::quantities>());
        s.executor().forall(ki, mdpolicy_zx) {
          auto [k, i] = ki;
          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto yhigh = bm[1][HIGH];
          if(yhigh == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m, k) = r(i, j - 1 - ghost_zone_size, k);
              v(i, j - 1 - m, k) = v(i, j - 1 - ghost_zone_size, k);
              p(i, j - 1 - m, k) = p(i, j - 1 - ghost_zone_size, k);
              sie(i, j - 1 - m, k) = sie(i, j - 1 - ghost_zone_size, k);

              ru(i, j - 1 - m, k) = ru(i, j - 1 - ghost_zone_size, k);
              rE(i, j - 1 - m, k) = rE(i, j - 1 - ghost_zone_size, k);
              Erad(i, j - 1 - m, k) = Erad(i, j - 1 - ghost_zone_size, k);
            }
          }
          else if(yhigh == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m, k) = r(i, j - 2 * ghost_zone_size + m, k);
              v(i, j - 1 - m, k).x = v(i, j - 2 * ghost_zone_size + m, k).x;
              v(i, j - 1 - m, k).y = -v(i, j - 2 * ghost_zone_size + m, k).y;
              v(i, j - 1 - m, k).z = v(i, j - 2 * ghost_zone_size + m, k).z;
              p(i, j - 1 - m, k) = p(i, j - 2 * ghost_zone_size + m, k);
              sie(i, j - 1 - m, k) = sie(i, j - 2 * ghost_zone_size + m, k);

              ru(i, j - 1 - m, k).x = ru(i, j - 2 * ghost_zone_size + m, k).x;
              ru(i, j - 1 - m, k).y = -ru(i, j - 2 * ghost_zone_size + m, k).y;
              ru(i, j - 1 - m, k).z = ru(i, j - 2 * ghost_zone_size + m, k).z;
              rE(i, j - 1 - m, k) = rE(i, j - 2 * ghost_zone_size + m, k);
              Erad(i, j - 1 - m, k) = Erad(i, j - 2 * ghost_zone_size + m, k);
            }
          }
        }; // forall
      } // if
    } // if
  }
  if constexpr(D == 3) {
    if(m.template is_low<ax::z>()) {
      auto mdpolicy_yx = get_mdiota_policy(r,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());
      s.executor().forall(ji, mdpolicy_yx) {
        auto [j, i] = ji;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto zlow = bm[2][LOW];
        if(zlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, m) = r(i, j, ghost_zone_size);
            v(i, j, m) = v(i, j, ghost_zone_size);
            p(i, j, m) = p(i, j, ghost_zone_size);
            sie(i, j, m) = sie(i, j, ghost_zone_size);

            ru(i, j, m) = ru(i, j, ghost_zone_size);
            rE(i, j, m) = rE(i, j, ghost_zone_size);
            Erad(i, j, m) = Erad(i, j, ghost_zone_size);
          }
        }
        else if(zlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, m) = r(i, j, 2 * ghost_zone_size - 1 - m);
            v(i, j, m).x = v(i, j, 2 * ghost_zone_size - 1 - m).x;
            v(i, j, m).y = v(i, j, 2 * ghost_zone_size - 1 - m).y;
            v(i, j, m).z = -v(i, j, 2 * ghost_zone_size - 1 - m).z;
            p(i, j, m) = p(i, j, 2 * ghost_zone_size - 1 - m);
            sie(i, j, m) = sie(i, j, 2 * ghost_zone_size - 1 - m);

            ru(i, j, m).x = ru(i, j, 2 * ghost_zone_size - 1 - m).x;
            ru(i, j, m).y = ru(i, j, 2 * ghost_zone_size - 1 - m).y;
            ru(i, j, m).z = -ru(i, j, 2 * ghost_zone_size - 1 - m).z;
            rE(i, j, m) = rE(i, j, 2 * ghost_zone_size - 1 - m);
            Erad(i, j, m) = Erad(i, j, 2 * ghost_zone_size - 1 - m);
          }
        }
      }; // forall
    }
    if(m.template is_high<ax::z>()) {
      const std::size_t k = m.template size<ax::z, dm::all>();
      auto mdpolicy_yx = get_mdiota_policy(r,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());
      s.executor().forall(ji, mdpolicy_yx) {
        auto [j, i] = ji;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto zhigh = bm[2][HIGH];
        if(zhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, k - 1 - m) = r(i, j, k - 1 - ghost_zone_size);
            v(i, j, k - 1 - m) = v(i, j, k - 1 - ghost_zone_size);
            p(i, j, k - 1 - m) = p(i, j, k - 1 - ghost_zone_size);
            sie(i, j, k - 1 - m) = sie(i, j, k - 1 - ghost_zone_size);

            ru(i, j, k - 1 - m) = ru(i, j, k - 1 - ghost_zone_size);
            rE(i, j, k - 1 - m) = rE(i, j, k - 1 - ghost_zone_size);
            Erad(i, j, k - 1 - m) = Erad(i, j, k - 1 - ghost_zone_size);
          }
        }
        else if(zhigh == bd::reflecting) {

          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, k - 1 - m) = r(i, j, k - 2 * ghost_zone_size + m);
            v(i, j, k - 1 - m).x = v(i, j, k - 2 * ghost_zone_size + m).x;
            v(i, j, k - 1 - m).y = v(i, j, k - 2 * ghost_zone_size + m).y;
            v(i, j, k - 1 - m).z = -v(i, j, k - 2 * ghost_zone_size + m).z;
            p(i, j, k - 1 - m) = p(i, j, k - 2 * ghost_zone_size + m);
            sie(i, j, k - 1 - m) = sie(i, j, k - 2 * ghost_zone_size + m);

            ru(i, j, k - 1 - m).x = ru(i, j, k - 2 * ghost_zone_size + m).x;
            ru(i, j, k - 1 - m).y = ru(i, j, k - 2 * ghost_zone_size + m).y;
            ru(i, j, k - 1 - m).z = -ru(i, j, k - 2 * ghost_zone_size + m).z;
            rE(i, j, k - 1 - m) = rE(i, j, k - 2 * ghost_zone_size + m);
            Erad(i, j, k - 1 - m) = Erad(i, j, k - 2 * ghost_zone_size + m);
          }
        }
      }; // forall
    } // if
  } // if
} // flow

// TODO: Undo the hardcoded boundary
template<std::size_t D>
void
apply_radiation_boundary(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<rw, na> r_a,
  flecsi::future<double> radiation_boundary_f) {
  using hard::tasks::util::get_mdiota_policy;
  auto r = m.template mdcolex<is::cells>(r_a);

  const size_t ghost_zone_size = m.ghost_zone_size();
  double radiation_boundary = radiation_boundary_f.get();

  if(m.template is_high<ax::x>()) {
    const std::size_t i = m.template size<ax::x, dm::all>();
    if constexpr(D == 1) {
      flecsi::util::iota_view policy{0, 1}; // default execution space
      s.executor().forall(x, policy) {
        (void)x;
        for(size_t m = 0; m < ghost_zone_size; ++m) {
          r(i - 1 - m) = radiation_boundary;
        }
      }; // forall
    }
    else if constexpr(D == 2) {
      s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
        for(size_t m = 0; m < ghost_zone_size; ++m) {
          r(i - 1 - m, j) = radiation_boundary;
        }
      }; // /forall
    }
    else /* D == 3 */ {
      auto mdpolicy_zy = get_mdiota_policy(r,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>());

      s.executor().forall(kj, mdpolicy_zy) {
        auto [k, j] = kj;
        for(size_t m = 0; m < ghost_zone_size; ++m) {
          r(i - 1 - m, j, k) = radiation_boundary;
        }
      }; // forall
    } // if
  } // if
} // flow

// TODO: Refactor. apply_boundary overload for just one field.
template<std::size_t D>
void
apply_boundary_single_field(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  typename single<typename mesh<D>::bmap>::template accessor<ro> bmap_a,
  field<double>::accessor<rw, na> r_a) {
  using hard::tasks::util::get_mdiota_policy;
  auto r = m.template mdcolex<is::cells>(r_a);

  const size_t ghost_zone_size = m.ghost_zone_size();

  if(m.template is_low<ax::x>()) {
    if constexpr(D == 1) {
      flecsi::util::iota_view policy{0, 1}; // default execution space
      s.executor().forall(x, policy) {
        (void)x; // remove compiler warning for unused

        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xlow = bm[0][LOW];
        if(xlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m) = r(ghost_zone_size);
          }
        }
        else if(xlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m) = r(2 * ghost_zone_size - 1 - m);
          }
        }
      };
    }
    else if constexpr(D == 2) {
      s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xlow = bm[0][LOW];
        if(xlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j) = r(ghost_zone_size, j);
          }
        }
        else if(xlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j) = r(2 * ghost_zone_size - 1 - m, j);
          }
        }
      }; // for
    }
    else /* D == 3 */ {
      auto mdpolicy_zy = get_mdiota_policy(r,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>());

      s.executor().forall(kj, mdpolicy_zy) {
        auto [k, j] = kj;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xlow = bm[0][LOW];
        if(xlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j, k) = r(ghost_zone_size, j, k);
          }
        }
        else if(xlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(m, j, k) = r(2 * ghost_zone_size - 1 - m, j, k);
          }
        }
      }; // forall
    } // if
  }
  if(m.template is_high<ax::x>()) {
    const std::size_t i = m.template size<ax::x, dm::all>();
    if constexpr(D == 1) {
      flecsi::util::iota_view policy{0, 1}; // default execution space
      s.executor().forall(x, policy) {
        (void)x; // remove compiler warning for unused

        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xhigh = bm[0][HIGH];
        if(xhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m) = r(i - 1 - ghost_zone_size);
          }
        }
        else if(xhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m) = r(i - 2 * ghost_zone_size + m);
          }
        }
      };
    }
    else if constexpr(D == 2) {
      s.executor().forall(j, (m.template cells<ax::y, dm::quantities>())) {
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xhigh = bm[0][HIGH];
        if(xhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j) = r(i - 1 - ghost_zone_size, j);
          }
        }
        else if(xhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j) = r(i - 2 * ghost_zone_size + m, j);
          }
        }
      }; // forall
    }
    else /* D == 3 */ {
      auto mdpolicy_kj = get_mdiota_policy(r,
        m.template cells<ax::z, dm::quantities>(),
        m.template cells<ax::y, dm::quantities>());

      s.executor().forall(kj, mdpolicy_kj) {
        auto [k, j] = kj;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto xhigh = bm[0][HIGH];
        if(xhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j, k) = r(i - 1 - ghost_zone_size, j, k);
          }
        }
        else if(xhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i - 1 - m, j, k) = r(i - 2 * ghost_zone_size + m, j, k);
          }
        }
      }; // forall
    } // if
  } // if

  if constexpr(D == 2 || D == 3) {
    if(m.template is_low<ax::y>()) {
      if constexpr(D == 2) {
        s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto ylow = bm[1][LOW];
          if(ylow == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m) = r(i, ghost_zone_size);
            }
          }
          else if(ylow == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m) = r(i, 2 * ghost_zone_size - 1 - m);
            }
          }
        }; // forall
      }
      else /* D == 3 */ {
        auto mdpolicy_zx = get_mdiota_policy(r,
          m.template cells<ax::z, dm::quantities>(),
          m.template cells<ax::x, dm::quantities>());
        s.executor().forall(ki, mdpolicy_zx) {
          auto [k, i] = ki;

          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto ylow = bm[1][LOW];
          if(ylow == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m, k) = r(i, ghost_zone_size, k);
            }
          }
          else if(ylow == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, m, k) = r(i, 2 * ghost_zone_size - 1 - m, k);
            }
          }
        }; // forall
      } // if
    }
    if(m.template is_high<ax::y>()) {
      const std::size_t j = m.template size<ax::y, dm::all>();
      if constexpr(D == 2) {
        s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto yhigh = bm[1][HIGH];
          if(yhigh == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m) = r(i, j - 1 - ghost_zone_size);
            }
          }
          else if(yhigh == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m) = r(i, j - 2 * ghost_zone_size + m);
            }
          }
        }; // forall
      }
      else /* D == 3 */ {
        auto mdpolicy_zx = get_mdiota_policy(r,
          m.template cells<ax::z, dm::quantities>(),
          m.template cells<ax::x, dm::quantities>());
        s.executor().forall(ki, mdpolicy_zx) {
          auto [k, i] = ki;
          const typename mesh<D>::bmap & bm = *bmap_a;
          const auto yhigh = bm[1][HIGH];
          if(yhigh == bd::flow) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m, k) = r(i, j - 1 - ghost_zone_size, k);
            }
          }
          else if(yhigh == bd::reflecting) {
            for(size_t m = 0; m < ghost_zone_size; ++m) {
              r(i, j - 1 - m, k) = r(i, j - 2 * ghost_zone_size + m, k);
            }
          }
        }; // forall
      } // if
    } // if
  }
  if constexpr(D == 3) {
    if(m.template is_low<ax::z>()) {
      auto mdpolicy_yx = get_mdiota_policy(r,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());
      s.executor().forall(ji, mdpolicy_yx) {
        auto [j, i] = ji;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto zlow = bm[2][LOW];
        if(zlow == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, m) = r(i, j, ghost_zone_size);
          }
        }
        else if(zlow == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, m) = r(i, j, 2 * ghost_zone_size - 1 - m);
          }
        }
      }; // forall
    }
    if(m.template is_high<ax::z>()) {
      const std::size_t k = m.template size<ax::z, dm::all>();
      auto mdpolicy_yx = get_mdiota_policy(r,
        m.template cells<ax::y, dm::quantities>(),
        m.template cells<ax::x, dm::quantities>());
      s.executor().forall(ji, mdpolicy_yx) {
        auto [j, i] = ji;
        const typename mesh<D>::bmap & bm = *bmap_a;
        const auto zhigh = bm[2][HIGH];
        if(zhigh == bd::flow) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, k - 1 - m) = r(i, j, k - 1 - ghost_zone_size);
          }
        }
        else if(zhigh == bd::reflecting) {
          for(size_t m = 0; m < ghost_zone_size; ++m) {
            r(i, j, k - 1 - m) = r(i, j, k - 2 * ghost_zone_size + m);
          }
        }
      }; // forall
    } // if
  } // if
} // flow

} // namespace hard::tasks
