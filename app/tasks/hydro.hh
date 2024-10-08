#ifndef HARD_TASKS_HYDRO_HH
#define HARD_TASKS_HYDRO_HH

#include "../numerical_algorithms/riemann_solvers.hh"
#include "../state.hh"
#include "hydro/cons2prim.hh"
#include "utils.hh"

#include <limits>

namespace hard::tasks::hydro {

template<std::size_t D>
double
update_dtmin(typename mesh<D>::template accessor<ro> m,
  flecsi::future<double> lmax_f) {

  double lmax = lmax_f.get();
  if constexpr(D == 1) {
    return m.template delta<ax::x>() / lmax;
  }
  else if constexpr(D == 2) {
    return std::min(
      m.template delta<ax::x>() / lmax, m.template delta<ax::y>() / lmax);
  }
  else {
    return std::min(m.template delta<ax::x>() / lmax,
      std::min(
        m.template delta<ax::y>() / lmax, m.template delta<ax::z>() / lmax));
  } // if
} // update_dtmin

#define DEBUG_PRINT(Ms, Mm, Mr, Mru, MrE)                                      \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Ms << std::endl;                                                     \
    for(auto j : m.template cells<ax::y, dm::quantities>()) {                  \
      for(auto i : m.template cells<ax::x, dm::predictor>()) {                 \
        ss << Mr(i, j, 2) << " ";                                              \
      }                                                                        \
      ss << std::endl;                                                         \
    }                                                                          \
    flog(info) << ss.str() << std::endl;                                       \
  }                                                                            \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Ms << std::endl;                                                     \
    for(auto j : m.template cells<ax::y, dm::quantities>()) {                  \
      for(auto i : m.template cells<ax::x, dm::predictor>()) {                 \
        ss << Mru(i, j, 2) << " ";                                             \
      }                                                                        \
      ss << std::endl;                                                         \
    }                                                                          \
    flog(info) << ss.str() << std::endl;                                       \
  }                                                                            \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Ms << std::endl;                                                     \
    for(auto j : m.template cells<ax::y, dm::quantities>()) {                  \
      for(auto i : m.template cells<ax::x, dm::predictor>()) {                 \
        ss << MrE(i, j, 2) << " ";                                             \
      }                                                                        \
      ss << std::endl;                                                         \
    }                                                                          \
    flog(info) << ss.str() << std::endl;                                       \
  }

template<ax::axis A, typename L, std::size_t D>
void
advance(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<rw, ro> r_a,
  typename field<vec<D>>::template accessor<rw, ro> ru_a,
  field<double>::accessor<rw, ro> rE_a,
  field<double>::accessor<rw, ro> Erad_a,
  typename field<vec<D>>::template accessor<rw, ro> u_a,
  field<double>::accessor<rw, ro> p_a,
  field<double>::accessor<wo, ro> q_a,
  typename field<vec<D>>::template accessor<wo, ro> qu_a,
  field<double>::accessor<wo, ro> qE_a,
  field<double>::accessor<wo, ro> qErad_a,
  field<double>::accessor<wo, ro> dr_ds_a,
  typename field<vec<D>>::template accessor<wo, ro> du_ds_a,
  field<double>::accessor<wo, ro> dp_ds_a,
  field<double>::accessor<wo, ro> dErad_ds_a,
  field<double>::accessor<wo, ro> rTail_a,
  typename field<vec<D>>::template accessor<wo, ro> ruTail_a,
  field<double>::accessor<wo, ro> rETail_a,
  typename field<vec<D>>::template accessor<wo, ro> uTail_a,
  field<double>::accessor<wo, ro> pTail_a,
  field<double>::accessor<wo, ro> EradTail_a,
  field<double>::accessor<wo, ro> rHead_a,
  typename field<vec<D>>::template accessor<wo, ro> ruHead_a,
  field<double>::accessor<wo, ro> rEHead_a,
  typename field<vec<D>>::template accessor<wo, ro> uHead_a,
  field<double>::accessor<wo, ro> pHead_a,
  field<double>::accessor<wo, ro> EradHead_a,
  field<double>::accessor<wo, ro> rF_a,
  typename field<vec<D>>::template accessor<wo, ro> ruF_a,
  field<double>::accessor<wo, ro> rEF_a,
  field<double>::accessor<wo, ro> EradF_a,
  single<double>::accessor<ro> gamma_a,
  double dt) {
  auto r = m.template mdcolex<is::cells>(r_a);
  auto ru = m.template mdcolex<is::cells>(ru_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);
  auto Erad = m.template mdcolex<is::cells>(Erad_a);
  auto u = m.template mdcolex<is::cells>(u_a);
  auto p = m.template mdcolex<is::cells>(p_a);
  auto q = m.template mdcolex<is::cells>(q_a);
  auto qu = m.template mdcolex<is::cells>(qu_a);
  auto qE = m.template mdcolex<is::cells>(qE_a);
  auto qErad = m.template mdcolex<is::cells>(qErad_a);
  auto dr_ds = m.template mdcolex<is::cells>(dr_ds_a);
  auto du_ds = m.template mdcolex<is::cells>(du_ds_a);
  auto dp_ds = m.template mdcolex<is::cells>(dp_ds_a);
  auto dErad_ds = m.template mdcolex<is::cells>(dErad_ds_a);
  auto rTail = m.template mdcolex<is::cells>(rTail_a);
  auto ruTail = m.template mdcolex<is::cells>(ruTail_a);
  auto rETail = m.template mdcolex<is::cells>(rETail_a);
  auto EradTail = m.template mdcolex<is::cells>(EradTail_a);
  auto uTail = m.template mdcolex<is::cells>(uTail_a);
  auto pTail = m.template mdcolex<is::cells>(pTail_a);
  auto rHead = m.template mdcolex<is::cells>(rHead_a);
  auto ruHead = m.template mdcolex<is::cells>(ruHead_a);
  auto rEHead = m.template mdcolex<is::cells>(rEHead_a);
  auto EradHead = m.template mdcolex<is::cells>(EradHead_a);
  auto uHead = m.template mdcolex<is::cells>(uHead_a);
  auto pHead = m.template mdcolex<is::cells>(pHead_a);
  auto rF = m.template mdcolex<is::cells>(rF_a);
  auto ruF = m.template mdcolex<is::cells>(ruF_a);
  auto rEF = m.template mdcolex<is::cells>(rEF_a);
  auto EradF = m.template mdcolex<is::cells>(EradF_a);
  auto const gamma = *gamma_a;

  if constexpr(D == 1) {
    /*------------------------------------------------------------------------*
      Compute slopes and reconstruct faces.
     *------------------------------------------------------------------------*/

    const double slope_factor{1.0 / m.template delta<A>()};
    const double extrap_factor{m.template delta<A>() / 2.0};
    const double p2rEFactor{1.0 / (gamma - 1.0)};
    const double courant{dt / (2.0 * m.template delta<A>())};
    const double gamma_minus_one = {gamma - 1.0};

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_slope_face_r_1d") {
      // Compute slope
      dr_ds(i) = slope_factor * L::limit(r(i - 1), r(i), r(i + 1));
      // Primitive quantities
      rTail(i) = r(i) + extrap_factor * dr_ds(i);
      rHead(i) = r(i) - extrap_factor * dr_ds(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_slope_face_u_1d") {
      // Compute slope
      du_ds(i).x = slope_factor * L::limit(u(i - 1).x, u(i).x, u(i + 1).x);
      // Primitive quantities
      uTail(i) = u(i) + extrap_factor * du_ds(i);
      uHead(i) = u(i) - extrap_factor * du_ds(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_face_ruTail_1d") {
      // Conserved quantities
      ruTail(i) = rTail(i) * uTail(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_face_ruHead_1d") {
      // Conserved quantities
      ruHead(i) = rHead(i) * uHead(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_slope_face_p_1d") {
      // Compute slope
      dp_ds(i) = slope_factor * L::limit(p(i - 1), p(i), p(i + 1));
      // Primitive quantities
      pTail(i) = p(i) + extrap_factor * dp_ds(i);
      pHead(i) = p(i) - extrap_factor * dp_ds(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_face_rETail_1d") {
      // Conserved quantities
      rETail(i) =
        p2rEFactor * pTail(i) + 0.5 * rTail(i) * utils::sqr(uTail(i).x);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::predictor>()), "adv_face_rEHead_1d") {
      // Conserved quantities
      rEHead(i) =
        p2rEFactor * pHead(i) + 0.5 * rHead(i) * utils::sqr(uHead(i).x);
    }; // for

    /*------------------------------------------------------------------------*
      Predictor step: approximate flux by averaging the governing flux
      functions at either side of the cell.
     *------------------------------------------------------------------------*/

    forall(i, (m.template cells<ax::x, dm::predictor>()), "adv_predictor_1d") {
      // Density
      q(i) = r(i) - courant * (ruTail(i).x - ruHead(i).x);

      // Momentum
      // ru^2 + p
      const vec<1> f_qu_Tail{ruTail(i).x * uTail(i).x + pTail(i)};
      const vec<1> f_qu_Head{ruHead(i).x * uHead(i).x + pHead(i)};
      qu(i) = ru(i) - courant * (f_qu_Tail - f_qu_Head);

      // Total Energy
      const double f_qE_Tail = (rETail(i) + pTail(i)) * uTail(i).x;
      const double f_qE_Head = (rEHead(i) + pHead(i)) * uHead(i).x;
      qE(i) = rE(i) - courant * (f_qE_Tail - f_qE_Head);
    }; // for

    forall(i, (m.template cells<ax::x, dm::predictor>()), "adv_primitives_1d") {
      // Primitives
      u(i) = qu(i) / q(i);
      p(i) = gamma_minus_one * (qE(i) - 0.5 * q(i) * utils::sqr(u(i).x));
    }; // for

      /*------------------------------------------------------------------------*
        Radiation slope face energy
       *------------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

    forall(i,
      (m.template cells<ax::x, dm::predictor>()),
      "adv_slope_face_energy_rad_1d") {
      dErad_ds(i) = slope_factor * L::limit(Erad(i - 1), Erad(i), Erad(i + 1));

      // Conserved quantities
      EradTail(i) = Erad(i) + extrap_factor * dErad_ds(i);
      EradHead(i) = Erad(i) - extrap_factor * dErad_ds(i);

      // Radiation energy
      qErad(i) = Erad(i) - courant * ((EradTail(i) * uTail(i).x) -
                                       (EradHead(i) * uHead(i).x));
    }; // for

    forall(i, (m.template cells<ax::x, dm::corrector>()), "adv_quant_rad_1d") {
      // Conserved quantities
      EradTail(i) = qErad(i - 1) + extrap_factor * dErad_ds(i - 1);
      EradHead(i) = qErad(i) - extrap_factor * dErad_ds(i);
    }; // for
#endif

    /*------------------------------------------------------------------------*
      Reconstruct faces from intermediates.
     *------------------------------------------------------------------------*/

    forall(i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_r_1d") {
      // Primitive quantities
      rTail(i) = q(i - 1) + extrap_factor * dr_ds(i - 1);
      rHead(i) = q(i) - extrap_factor * dr_ds(i);
    }; // for

    forall(i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_u_1d") {
      // Primitive quantities
      uTail(i) = u(i - 1) + extrap_factor * du_ds(i - 1);
      uHead(i) = u(i) - extrap_factor * du_ds(i);
    }; // for

    forall(i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_p_1d") {
      // Primitive quantities
      pTail(i) = p(i - 1) + extrap_factor * dp_ds(i - 1);
      pHead(i) = p(i) - extrap_factor * dp_ds(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_ruTail_1d") {
      // Conserved quantities
      ruTail(i) = rTail(i) * uTail(i);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_ruHead_1d") {
      // Conserved quantities
      ruHead(i).x = rHead(i) * uHead(i).x;
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_rETail_1d") {
      // Conserved quantities
      rETail(i) =
        p2rEFactor * pTail(i) + 0.5 * rTail(i) * utils::sqr(uTail(i).x);
    }; // for

    forall(
      i, (m.template cells<ax::x, dm::corrector>()), "adv_re_face_rEHead_1d") {
      // Conserved quantities
      rEHead(i) =
        p2rEFactor * pHead(i) + 0.5 * rHead(i) * utils::sqr(uHead(i).x);
    }; // for

    /*------------------------------------------------------------------------*
      Compute corrector setp: Riemann solver.
     *------------------------------------------------------------------------*/

    forall(i, (m.template cells<ax::x, dm::corrector>()), "adv_riemann_1d") {

      // Update min/max eigenvalues on tail face.
      const double cT = std::sqrt(gamma * pTail(i) / rTail(i));
      const double LminT = uTail(i).x - cT;
      const double LmaxT = uTail(i).x + cT;

      // Update min/max eigenvalues on head face.
      const double cH = std::sqrt(gamma * pHead(i) / rHead(i));
      const double LminH = uHead(i).x - cH;
      const double LmaxH = uHead(i).x + cH;

      // Fluxes from left and right state, for the Riemann solver
      const double f_r_T{ruTail(i).x};
      const double f_r_H{ruHead(i).x};
      const vec<1> f_ru_T{ruTail(i).x * uTail(i).x + pTail(i)};
      const vec<1> f_ru_H{ruHead(i).x * uHead(i).x + pHead(i)};
      const double f_rE_T{(rETail(i) + pTail(i)) * uTail(i).x};
      const double f_rE_H{(rEHead(i) + pHead(i)) * uHead(i).x};

      const auto hll_fluxes = numerical_algorithms::hll_hydro(rTail(i),
        rHead(i),
        ruTail(i),
        ruHead(i),
        rETail(i),
        rEHead(i),
        f_r_T,
        f_r_H,
        f_ru_T,
        f_ru_H,
        f_rE_T,
        f_rE_H,
        LminT,
        LmaxT,
        LminH,
        LmaxH);

      rF(i) = get<0>(hll_fluxes);
      ruF(i) = get<1>(hll_fluxes);
      rEF(i) = get<2>(hll_fluxes);
    }; // for

    /*------------------------------------------------------------------------*
      Update averages (in-place).
     *------------------------------------------------------------------------*/

    const double update_factor = dt / m.template delta<A>();
    forall(i, (m.template cells<ax::x, dm::quantities>()), "adv_upd_r_1d") {
      r(i) = r(i) + update_factor * (rF(i) - rF(i + 1));
    }; // for

    forall(i, (m.template cells<ax::x, dm::quantities>()), "adv_upd_ru_1d") {
      ru(i) = ru(i) + update_factor * (ruF(i) - ruF(i + 1));
    }; // for

    forall(i, (m.template cells<ax::x, dm::quantities>()), "adv_upd_rE_1d") {
      rE(i) = rE(i) + update_factor * (rEF(i) - rEF(i + 1));
    }; // for

      /*------------------------------------------------------------------------*
        Radiation Update
       *------------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

    forall(i, (m.template cells<ax::x, dm::corrector>()), "adv_advect_rad_1d") {
      EradF(i) = numerical_algorithms::advect_Erad(
        EradTail(i), EradHead(i), uTail(i).x, uHead(i).x);
    }; // for

    forall(i, (m.template cells<ax::x, dm::quantities>()), "adv_upd_Erad_1d") {
      Erad(i) = Erad(i) + update_factor * (EradF(i) - EradF(i + 1));
    }; // for

#endif
  }
  else if constexpr(D == 2) {
    if constexpr(A == ax::x) {
      const double slope_factor{1.0 / m.template delta<A>()};
      const double extrap_factor{m.template delta<A>() / 2.0};
      const double p2rEFactor{1.0 / (gamma - 1.0)};
      const double courant{dt / (2.0 * m.template delta<A>())};
      const double update_factor = dt / m.template delta<A>();

      forall(j,
        (m.template cells<ax::y, dm::quantities>()),
        "adv_slope_face_energy_2dx") {
        for(auto i : m.template cells<ax::x, dm::predictor>()) {
          /*----------------------------------------------------------------------*
           Compute slopes and reconstruct faces
          *----------------------------------------------------------------------*/
          // Compute slope
          dr_ds(i, j) =
            slope_factor * L::limit(r(i - 1, j), r(i, j), r(i + 1, j));
          // Primitive quantity face
          rTail(i, j) = r(i, j) + extrap_factor * dr_ds(i, j);
          rHead(i, j) = r(i, j) - extrap_factor * dr_ds(i, j);

          // Compute slope
          du_ds(i, j).x =
            slope_factor * L::limit(u(i - 1, j).x, u(i, j).x, u(i + 1, j).x);
          du_ds(i, j).y =
            slope_factor * L::limit(u(i - 1, j).y, u(i, j).y, u(i + 1, j).y);
          // Primitive quantity face
          uTail(i, j) = u(i, j) + extrap_factor * du_ds(i, j);
          uHead(i, j) = u(i, j) - extrap_factor * du_ds(i, j);

          // Conserved quantities face
          ruTail(i, j) = rTail(i, j) * uTail(i, j);
          ruHead(i, j) = rHead(i, j) * uHead(i, j);

          // Compute slope
          dp_ds(i, j) =
            slope_factor * L::limit(p(i - 1, j), p(i, j), p(i + 1, j));
          // Primitive quantity face
          pTail(i, j) = p(i, j) + extrap_factor * dp_ds(i, j);
          pHead(i, j) = p(i, j) - extrap_factor * dp_ds(i, j);

          // Conserved quantities face
          rETail(i, j) = p2rEFactor * pTail(i, j) +
                         0.5 * rTail(i, j) * uTail(i, j).norm_squared();
          rEHead(i, j) = p2rEFactor * pHead(i, j) +
                         0.5 * rHead(i, j) * uHead(i, j).norm_squared();

          /*----------------------------------------------------------------------*
           Predictor step: approximate flux by averaging the governing flux
           functions at either side of the cell.
          *----------------------------------------------------------------------*/
          // Density
          q(i, j) = r(i, j) - courant * (ruTail(i, j).x - ruHead(i, j).x);

          // Momentum
          // ru^2 + p
          // ruv
          const vec<2> f_qu_Tail{ruTail(i, j).x * uTail(i, j).x + pTail(i, j),
            rTail(i, j) * uTail(i, j).x * uTail(i, j).y};
          const vec<2> f_qu_Head{ruHead(i, j).x * uHead(i, j).x + pHead(i, j),
            rHead(i, j) * uHead(i, j).x * uHead(i, j).y};
          qu(i, j) = ru(i, j) - courant * (f_qu_Tail - f_qu_Head);

          // Total Energy
          const double f_qE_Tail = (rETail(i, j) + pTail(i, j)) * uTail(i, j).x;
          const double f_qE_Head = (rEHead(i, j) + pHead(i, j)) * uHead(i, j).x;
          qE(i, j) = rE(i, j) - courant * (f_qE_Tail - f_qE_Head);
        }
      };

      forall(j,
        (m.template cells<ax::y, dm::quantities>()),
        "adv_update_primitives_2dx") {
        for(auto i : m.template cells<ax::x, dm::predictor>()) {
          // Primitives
          u(i, j) = qu(i, j) / q(i, j);
          p(i, j) =
            (gamma - 1.0) * (qE(i, j) - 0.5 * q(i, j) * u(i, j).norm_squared());
        }
      };

        /*----------------------------------------------------------------------*
         Raddiation slope, faces, energy
        *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(j,
        (m.template cells<ax::y, dm::quantities>()),
        "adv_slope_face_energy_rad_2dx") {
        for(auto i : m.template cells<ax::x, dm::predictor>()) {
          // Compute slopes
          dErad_ds(i, j) =
            slope_factor * L::limit(Erad(i - 1, j), Erad(i, j), Erad(i + 1, j));
          // Conserved quantities
          EradTail(i, j) = Erad(i, j) + extrap_factor * dErad_ds(i, j);
          EradHead(i, j) = Erad(i, j) - extrap_factor * dErad_ds(i, j);
          // Radiation energy
          qErad(i, j) =
            Erad(i, j) - courant * ((EradTail(i, j) * uTail(i, j).x) -
                                     (EradHead(i, j) * uHead(i, j).x));
        } // for
      };
#endif

      forall(j,
        (m.template cells<ax::y, dm::quantities>()),
        "adv_face_riemann_2dx") {
        for(auto i : m.template cells<ax::x, dm::corrector>()) {
          /*----------------------------------------------------------------------*
           Reconstruct faces from intermediates.
          *----------------------------------------------------------------------*/
          rTail(i, j) = q(i - 1, j) + extrap_factor * dr_ds(i - 1, j);
          rHead(i, j) = q(i, j) - extrap_factor * dr_ds(i, j);
          uTail(i, j) = u(i - 1, j) + extrap_factor * du_ds(i - 1, j);
          uHead(i, j) = u(i, j) - extrap_factor * du_ds(i, j);
          // Conserved quantities
          ruTail(i, j) = rTail(i, j) * uTail(i, j);
          ruHead(i, j) = rHead(i, j) * uHead(i, j);

          pTail(i, j) = p(i - 1, j) + extrap_factor * dp_ds(i - 1, j);
          pHead(i, j) = p(i, j) - extrap_factor * dp_ds(i, j);
          // Conserved quantities
          rETail(i, j) = p2rEFactor * pTail(i, j) +
                         0.5 * rTail(i, j) * uTail(i, j).norm_squared();
          rEHead(i, j) = p2rEFactor * pHead(i, j) +
                         0.5 * rHead(i, j) * uHead(i, j).norm_squared();
          /*----------------------------------------------------------------------*
          Compute corrector setp: Riemann solver.
          *----------------------------------------------------------------------*/
          // Update min/max eigenvalues on tail face.
          const double cT = std::sqrt(gamma * pTail(i, j) / rTail(i, j));
          const double LminT = uTail(i, j).x - cT;
          const double LmaxT = uTail(i, j).x + cT;

          // Update min/max eigenvalues on head face.
          const double cH = std::sqrt(gamma * pHead(i, j) / rHead(i, j));
          const double LminH = uHead(i, j).x - cH;
          const double LmaxH = uHead(i, j).x + cH;

          // Fluxes from left and right state, for the Riemann solver
          const double f_r_T{ruTail(i, j).x};
          const double f_r_H{ruHead(i, j).x};
          const vec<2> f_ru_T{ruTail(i, j).x * uTail(i, j).x + pTail(i, j),
            ruTail(i, j).x * uTail(i, j).y};
          const vec<2> f_ru_H{ruHead(i, j).x * uHead(i, j).x + pHead(i, j),
            ruHead(i, j).x * uHead(i, j).y};
          const double f_rE_T{(rETail(i, j) + pTail(i, j)) * uTail(i, j).x};
          const double f_rE_H{(rEHead(i, j) + pHead(i, j)) * uHead(i, j).x};

          const auto hll_fluxes = numerical_algorithms::hll_hydro(rTail(i, j),
            rHead(i, j),
            ruTail(i, j),
            ruHead(i, j),
            rETail(i, j),
            rEHead(i, j),
            f_r_T,
            f_r_H,
            f_ru_T,
            f_ru_H,
            f_rE_T,
            f_rE_H,
            LminT,
            LmaxT,
            LminH,
            LmaxH);

          rF(i, j) = get<0>(hll_fluxes);
          ruF(i, j) = get<1>(hll_fluxes);
          rEF(i, j) = get<2>(hll_fluxes);
        } // for
      };

      /*----------------------------------------------------------------------*
        Update averages (in-place).
      *----------------------------------------------------------------------*/

      forall(j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_r_2dx") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          r(i, j) = r(i, j) + update_factor * (rF(i, j) - rF(i + 1, j));
        } // for
      };

      forall(j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_ru_2dx") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          ru(i, j) = ru(i, j) + update_factor * (ruF(i, j) - ruF(i + 1, j));
        } // for
      };

      forall(j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_rE_2dx") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          rE(i, j) = rE(i, j) + update_factor * (rEF(i, j) - rEF(i + 1, j));
        } // for
      };

        /*----------------------------------------------------------------------*
          Radiation update
        *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(
        j, (m.template cells<ax::y, dm::quantities>()), "adv_advect_Erad_2dx") {
        for(auto i : m.template cells<ax::x, dm::corrector>()) {
          // Conserved quantities
          EradTail(i, j) = qErad(i - 1, j) + extrap_factor * dErad_ds(i - 1, j);
          EradHead(i, j) = qErad(i, j) - extrap_factor * dErad_ds(i, j);

          EradF(i, j) = numerical_algorithms::advect_Erad(
            EradTail(i, j), EradHead(i, j), uTail(i, j).x, uHead(i, j).x);
        } // for
      };

      forall(
        j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_Erad_2dx") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          Erad(i, j) =
            Erad(i, j) + update_factor * (EradF(i, j) - EradF(i + 1, j));
        } // for
      };

#endif
    }
    else /* A == ax::y */ {
      const double slope_factor{1.0 / m.template delta<A>()};
      const double extrap_factor{m.template delta<A>() / 2.0};
      const double p2rEFactor{1.0 / (gamma - 1.0)};
      const double courant{dt / (2.0 * m.template delta<A>())};
      const double update_factor = dt / m.template delta<A>();

      forall(j,
        (m.template cells<ax::y, dm::predictor>()),
        "adv_slope_face_energy_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          /*----------------------------------------------------------------------*
            Compute slopes and reconstruct faces
          *----------------------------------------------------------------------*/
          // Compute slope
          dr_ds(i, j) =
            slope_factor * L::limit(r(i, j - 1), r(i, j), r(i, j + 1));
          // Primitive quantities
          rTail(i, j) = r(i, j) + extrap_factor * dr_ds(i, j);
          rHead(i, j) = r(i, j) - extrap_factor * dr_ds(i, j);

          // Compute slope
          du_ds(i, j).x =
            slope_factor * L::limit(u(i, j - 1).x, u(i, j).x, u(i, j + 1).x);
          du_ds(i, j).y =
            slope_factor * L::limit(u(i, j - 1).y, u(i, j).y, u(i, j + 1).y);
          // Primitive quantities
          uTail(i, j) = u(i, j) + extrap_factor * du_ds(i, j);
          uHead(i, j) = u(i, j) - extrap_factor * du_ds(i, j);

          // Conserved quantities
          ruTail(i, j) = rTail(i, j) * uTail(i, j);
          ruHead(i, j) = rHead(i, j) * uHead(i, j);

          // Compute slope
          dp_ds(i, j) =
            slope_factor * L::limit(p(i, j - 1), p(i, j), p(i, j + 1));
          // Primitive quantities
          pTail(i, j) = p(i, j) + extrap_factor * dp_ds(i, j);
          pHead(i, j) = p(i, j) - extrap_factor * dp_ds(i, j);

          // Conserved quantities
          rETail(i, j) = p2rEFactor * pTail(i, j) +
                         0.5 * rTail(i, j) * uTail(i, j).norm_squared();
          rEHead(i, j) = p2rEFactor * pHead(i, j) +
                         0.5 * rHead(i, j) * uHead(i, j).norm_squared();

          /*----------------------------------------------------------------------*
            Predictor step: approximate flux by averaging the governing flux
            functions at either side of the cell.
          *----------------------------------------------------------------------*/
          // Density
          q(i, j) = r(i, j) - courant * (ruTail(i, j).y - ruHead(i, j).y);

          // Momentum
          // ruv
          // rv^2+p
          // rvw
          const vec<2> f_qu_Tail{rTail(i, j) * uTail(i, j).x * uTail(i, j).y,
            ruTail(i, j).y * uTail(i, j).y + pTail(i, j)};
          const vec<2> f_qu_Head{rHead(i, j) * uHead(i, j).x * uHead(i, j).y,
            ruHead(i, j).y * uHead(i, j).y + pHead(i, j)};
          qu(i, j) = ru(i, j) - courant * (f_qu_Tail - f_qu_Head);

          // Total Energy
          const double f_qE_Tail = (rETail(i, j) + pTail(i, j)) * uTail(i, j).y;
          const double f_qE_Head = (rEHead(i, j) + pHead(i, j)) * uHead(i, j).y;
          qE(i, j) = rE(i, j) - courant * (f_qE_Tail - f_qE_Head);
        } // for
      }; // for

      forall(j,
        (m.template cells<ax::y, dm::predictor>()),
        "adv_update_primitives_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          // Primitives
          u(i, j) = qu(i, j) / q(i, j);
          p(i, j) =
            (gamma - 1.0) * (qE(i, j) - 0.5 * q(i, j) * u(i, j).norm_squared());
        } // for
      }; // for

        /*----------------------------------------------------------------------*
         Raddiation slope, faces, energy
        *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(j,
        (m.template cells<ax::y, dm::predictor>()),
        "adv_slope_face_energy_rad_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          // Compute Slope
          dErad_ds(i, j) =
            slope_factor * L::limit(Erad(i, j - 1), Erad(i, j), Erad(i, j + 1));
          // Conserved quantities
          EradTail(i, j) = Erad(i, j) + extrap_factor * dErad_ds(i, j);
          EradHead(i, j) = Erad(i, j) - extrap_factor * dErad_ds(i, j);
          // Radiation energy
          qErad(i, j) =
            Erad(i, j) - courant * ((EradTail(i, j) * uTail(i, j).y) -
                                     (EradHead(i, j) * uHead(i, j).y));
        } // for
      }; // for
#endif

      forall(
        j, (m.template cells<ax::y, dm::corrector>()), "adv_face_riemann_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          /*----------------------------------------------------------------------*
            Reconstruct faces from intermediates.
          *----------------------------------------------------------------------*/ // Primitive quantities
          rTail(i, j) = q(i, j - 1) + extrap_factor * dr_ds(i, j - 1);
          rHead(i, j) = q(i, j) - extrap_factor * dr_ds(i, j);
          uTail(i, j) = u(i, j - 1) + extrap_factor * du_ds(i, j - 1);
          uHead(i, j) = u(i, j) - extrap_factor * du_ds(i, j);
          // Conserved quantities
          ruTail(i, j) = rTail(i, j) * uTail(i, j);
          ruHead(i, j) = rHead(i, j) * uHead(i, j);

          pTail(i, j) = p(i, j - 1) + extrap_factor * dp_ds(i, j - 1);
          pHead(i, j) = p(i, j) - extrap_factor * dp_ds(i, j);
          // Conserved quantities
          rETail(i, j) = p2rEFactor * pTail(i, j) +
                         0.5 * rTail(i, j) * uTail(i, j).norm_squared();
          rEHead(i, j) = p2rEFactor * pHead(i, j) +
                         0.5 * rHead(i, j) * uHead(i, j).norm_squared();

          /*----------------------------------------------------------------------*
            Compute corrector setp: Riemann solver.
          *----------------------------------------------------------------------*/
          // Update min/max eigenvalues on tail face.
          const double cT = std::sqrt(gamma * pTail(i, j) / rTail(i, j));
          const double LminT = uTail(i, j).y - cT;
          const double LmaxT = uTail(i, j).y + cT;

          // Update min/max eigenvalues on head face.
          const double cH = std::sqrt(gamma * pHead(i, j) / rHead(i, j));
          const double LminH = uHead(i, j).y - cH;
          const double LmaxH = uHead(i, j).y + cH;

          // Fluxes from left and right state, for the Riemann solver
          const double f_r_T{ruTail(i, j).y};
          const double f_r_H{ruHead(i, j).y};
          const vec<2> f_ru_T{rTail(i, j) * uTail(i, j).x * uTail(i, j).y,
            ruTail(i, j).y * uTail(i, j).y + pTail(i, j)};
          const vec<2> f_ru_H{rHead(i, j) * uHead(i, j).x * uHead(i, j).y,
            ruHead(i, j).y * uHead(i, j).y + pHead(i, j)};
          const double f_rE_T{(rETail(i, j) + pTail(i, j)) * uTail(i, j).y};
          const double f_rE_H{(rEHead(i, j) + pHead(i, j)) * uHead(i, j).y};

          const auto hll_fluxes = numerical_algorithms::hll_hydro(rTail(i, j),
            rHead(i, j),
            ruTail(i, j),
            ruHead(i, j),
            rETail(i, j),
            rEHead(i, j),
            f_r_T,
            f_r_H,
            f_ru_T,
            f_ru_H,
            f_rE_T,
            f_rE_H,
            LminT,
            LmaxT,
            LminH,
            LmaxH);

          rF(i, j) = get<0>(hll_fluxes);
          ruF(i, j) = get<1>(hll_fluxes);
          rEF(i, j) = get<2>(hll_fluxes);
        } // for
      }; // for

      /*----------------------------------------------------------------------*
        Update averages (in-place).
       *----------------------------------------------------------------------*/

      forall(j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_r_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          r(i, j) = r(i, j) + update_factor * (rF(i, j) - rF(i, j + 1));
        } // for
      }; // for

      forall(j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_ru_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          ru(i, j) = ru(i, j) + update_factor * (ruF(i, j) - ruF(i, j + 1));
        } // for
      }; // for

      forall(j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_rE_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          rE(i, j) = rE(i, j) + update_factor * (rEF(i, j) - rEF(i, j + 1));
        } // for
      }; // for

        /*----------------------------------------------------------------------*
          Radiation update
        *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(
        j, (m.template cells<ax::y, dm::corrector>()), "adv_advect_Erad_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          // Conserved quantities
          EradTail(i, j) = qErad(i, j - 1) + extrap_factor * dErad_ds(i, j - 1);
          EradHead(i, j) = qErad(i, j) - extrap_factor * dErad_ds(i, j);

          EradF(i, j) = numerical_algorithms::advect_Erad(
            EradTail(i, j), EradHead(i, j), uTail(i, j).y, uHead(i, j).y);
        } // for
      }; // for

      forall(
        j, (m.template cells<ax::y, dm::quantities>()), "adv_upd_Erad_2dy") {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          Erad(i, j) =
            Erad(i, j) - update_factor * (EradF(i, j + 1) - EradF(i, j));
        } // for
      }; // for

#endif
    } // if
  }
  else /* D == 3 */ {
    if constexpr(A == ax::x) {

      const double slope_factor{1.0 / m.template delta<A>()};
      const double extrap_factor{m.template delta<A>() / 2.0};
      const double p2rEFactor{1.0 / (gamma - 1.0)};
      const double courant{dt / (2.0 * m.template delta<A>())};
      const double update_factor = dt / m.template delta<A>();

      forall(k,
        (m.template cells<ax::z, dm::quantities>()),
        "adv_slope_face_energy_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::predictor>()) {
            /*----------------------------------------------------------------------*
              Compute slopes, reconstruct faces.
            *----------------------------------------------------------------------*/
            // Slope
            dr_ds(i, j, k) =
              slope_factor *
              L::limit(r(i - 1, j, k), r(i, j, k), r(i + 1, j, k));
            // Primitive quantities
            rTail(i, j, k) = r(i, j, k) + extrap_factor * dr_ds(i, j, k);
            rHead(i, j, k) = r(i, j, k) - extrap_factor * dr_ds(i, j, k);

            // Slope
            du_ds(i, j, k).x =
              slope_factor *
              L::limit(u(i - 1, j, k).x, u(i, j, k).x, u(i + 1, j, k).x);
            du_ds(i, j, k).y =
              slope_factor *
              L::limit(u(i - 1, j, k).y, u(i, j, k).y, u(i + 1, j, k).y);
            du_ds(i, j, k).z =
              slope_factor *
              L::limit(u(i - 1, j, k).z, u(i, j, k).z, u(i + 1, j, k).z);
            // Primitive quantities
            uTail(i, j, k) = u(i, j, k) + extrap_factor * du_ds(i, j, k);
            uHead(i, j, k) = u(i, j, k) - extrap_factor * du_ds(i, j, k);

            // Conserved quantities
            ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
            ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);

            // Slope
            dp_ds(i, j, k) =
              slope_factor *
              L::limit(p(i - 1, j, k), p(i, j, k), p(i + 1, j, k));
            // Primitive quantities
            pTail(i, j, k) = p(i, j, k) + extrap_factor * dp_ds(i, j, k);
            pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

            // Conserved quantities
            rETail(i, j, k) =
              p2rEFactor * pTail(i, j, k) +
              0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
            rEHead(i, j, k) =
              p2rEFactor * pHead(i, j, k) +
              0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();

            /*----------------------------------------------------------------------*
             Predictor step: approximate flux by averaging the governing flux
             functions at either side of the cell.
            *----------------------------------------------------------------------*/
            // Density
            q(i, j, k) =
              r(i, j, k) - courant * (ruTail(i, j, k).x - ruHead(i, j, k).x);

            // Momentum
            // ru^2 + p
            // ruv
            // ruw
            const vec<3> f_qu_Tail{
              ruTail(i, j, k).x * uTail(i, j, k).x + pTail(i, j, k),
              rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).y,
              rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).z};
            const vec<3> f_qu_Head{
              ruHead(i, j, k).x * uHead(i, j, k).x + pHead(i, j, k),
              rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).y,
              rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).z};
            qu(i, j, k) = ru(i, j, k) - courant * (f_qu_Tail - f_qu_Head);

            // Total Energy
            const double f_qE_Tail =
              (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).x;
            const double f_qE_Head =
              (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).x;
            qE(i, j, k) = rE(i, j, k) - courant * (f_qE_Tail - f_qE_Head);
          } // for
        }
      };

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_primitives_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::predictor>()) {
            // Primitives
            u(i, j, k) = qu(i, j, k) / q(i, j, k);
            p(i, j, k) =
              (gamma - 1.0) *
              (qE(i, j, k) - 0.5 * q(i, j, k) * u(i, j, k).norm_squared());
          } // for
        }
      };

        /*----------------------------------------------------------------------*
          Radiation slope, face, and energy.
         *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(k,
        (m.template cells<ax::z, dm::quantities>()),
        "adv_slope_face_energy_rad_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::predictor>()) {
            // Slope
            dErad_ds(i, j, k) =
              slope_factor *
              L::limit(Erad(i - 1, j, k), Erad(i, j, k), Erad(i + 1, j, k));
            // Conserved quantities
            EradTail(i, j, k) =
              Erad(i, j, k) + extrap_factor * dErad_ds(i, j, k);
            EradHead(i, j, k) =
              Erad(i, j, k) - extrap_factor * dErad_ds(i, j, k);
            // Radiation Energy
            qErad(i, j, k) =
              Erad(i, j, k) -
              courant * ((EradTail(i, j, k) * uTail(i, j, k).x) -
                          (EradHead(i, j, k) * uHead(i, j, k).x));
          } // for
        }
      };

#endif

      forall(k,
        (m.template cells<ax::z, dm::quantities>()),
        "adv_re_face_riemann_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::corrector>()) {
            /*----------------------------------------------------------------------*
              Reconstruct faces from intermediates.
            *----------------------------------------------------------------------*/
            // Primitive quantities
            rTail(i, j, k) =
              q(i - 1, j, k) + extrap_factor * dr_ds(i - 1, j, k);
            rHead(i, j, k) = q(i, j, k) - extrap_factor * dr_ds(i, j, k);
            uTail(i, j, k) =
              u(i - 1, j, k) + extrap_factor * du_ds(i - 1, j, k);
            uHead(i, j, k) = u(i, j, k) - extrap_factor * du_ds(i, j, k);

            // Primitive quantities
            pTail(i, j, k) =
              p(i - 1, j, k) + extrap_factor * dp_ds(i - 1, j, k);
            pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

            // Conserved quantities
            ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
            ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);
            rETail(i, j, k) =
              p2rEFactor * pTail(i, j, k) +
              0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
            rEHead(i, j, k) =
              p2rEFactor * pHead(i, j, k) +
              0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();
            EradTail(i, j, k) =
              qErad(i - 1, j, k) + extrap_factor * dErad_ds(i - 1, j, k);
            EradHead(i, j, k) =
              qErad(i, j, k) - extrap_factor * dErad_ds(i, j, k);

            /*----------------------------------------------------------------------*
              Compute corrector setp: Riemann solver.
            *----------------------------------------------------------------------*/
            // Update min/max eigenvalues on tail face.
            const double cT =
              std::sqrt(gamma * pTail(i, j, k) / rTail(i, j, k));
            const double LminT = uTail(i, j, k).x - cT;
            const double LmaxT = uTail(i, j, k).x + cT;

            // Update min/max eigenvalues on head face.
            const double cH =
              std::sqrt(gamma * pHead(i, j, k) / rHead(i, j, k));
            const double LminH = uHead(i, j, k).x - cH;
            const double LmaxH = uHead(i, j, k).x + cH;

            // Fluxes from left and right state, for the Riemann solver
            const double f_r_T{ruTail(i, j, k).x};
            const double f_r_H{ruHead(i, j, k).x};
            const vec<3> f_ru_T{
              ruTail(i, j, k).x * uTail(i, j, k).x + pTail(i, j, k),
              ruTail(i, j, k).x * uTail(i, j, k).y,
              ruTail(i, j, k).x * uTail(i, j, k).z};
            const vec<3> f_ru_H{
              ruHead(i, j, k).x * uHead(i, j, k).x + pHead(i, j, k),
              ruHead(i, j, k).x * uHead(i, j, k).y,
              ruHead(i, j, k).x * uHead(i, j, k).z};
            const double f_rE_T{
              (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).x};
            const double f_rE_H{
              (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).x};

            const auto hll_fluxes =
              numerical_algorithms::hll_hydro(rTail(i, j, k),
                rHead(i, j, k),
                ruTail(i, j, k),
                ruHead(i, j, k),
                rETail(i, j, k),
                rEHead(i, j, k),
                f_r_T,
                f_r_H,
                f_ru_T,
                f_ru_H,
                f_rE_T,
                f_rE_H,
                LminT,
                LmaxT,
                LminH,
                LmaxH);

            rF(i, j, k) = get<0>(hll_fluxes);
            ruF(i, j, k) = get<1>(hll_fluxes);
            rEF(i, j, k) = get<2>(hll_fluxes);
            EradF(i, j, k) =
              numerical_algorithms::advect_Erad(EradTail(i, j, k),
                EradHead(i, j, k),
                uTail(i, j, k).x,
                uHead(i, j, k).x);

          } // for
        }
      };

      /*----------------------------------------------------------------------*
        Update averages (in-place).
       *----------------------------------------------------------------------*/

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_r_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            r(i, j, k) =
              r(i, j, k) + update_factor * (rF(i, j, k) - rF(i + 1, j, k));
          } // for
        }
      };

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_ru_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            ru(i, j, k) =
              ru(i, j, k) + update_factor * (ruF(i, j, k) - ruF(i + 1, j, k));
          } // for
        }
      };

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_rE_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            rE(i, j, k) =
              rE(i, j, k) + update_factor * (rEF(i, j, k) - rEF(i + 1, j, k));
          } // for
        }
      };

        /*----------------------------------------------------------------------*
          Radiation update (in-place).
         *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_re_face_rad_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::corrector>()) {
            // Conserved quantities
            EradTail(i, j, k) =
              qErad(i - 1, j, k) + extrap_factor * dErad_ds(i - 1, j, k);
            EradHead(i, j, k) =
              qErad(i, j, k) - extrap_factor * dErad_ds(i, j, k);
          } // for
        }
      };

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_Erad_3dx") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            Erad(i, j, k) =
              Erad(i, j, k) +
              update_factor * (EradF(i, j, k) - EradF(i + 1, j, k));
          } // for
        }
      };

#endif
    }
    else if constexpr(A == ax::y) {
      const double slope_factor{1.0 / m.template delta<A>()};
      const double extrap_factor{m.template delta<A>() / 2.0};
      const double p2rEFactor{1.0 / (gamma - 1.0)};
      const double courant{dt / (2.0 * m.template delta<A>())};
      const double update_factor = dt / m.template delta<A>();

      forall(k,
        (m.template cells<ax::z, dm::quantities>()),
        "adv_slope_face_energy_3dy") {
        for(auto j : m.template cells<ax::y, dm::predictor>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            /*----------------------------------------------------------------------*
              Compute slopes and reconstruct faces.
            *----------------------------------------------------------------------*/
            // Compute slope
            dr_ds(i, j, k) =
              slope_factor *
              L::limit(r(i, j - 1, k), r(i, j, k), r(i, j + 1, k));
            // Primitive quantities
            rTail(i, j, k) = r(i, j, k) + extrap_factor * dr_ds(i, j, k);
            rHead(i, j, k) = r(i, j, k) - extrap_factor * dr_ds(i, j, k);

            // Compute slope
            du_ds(i, j, k).x =
              slope_factor *
              L::limit(u(i, j - 1, k).x, u(i, j, k).x, u(i, j + 1, k).x);
            du_ds(i, j, k).y =
              slope_factor *
              L::limit(u(i, j - 1, k).y, u(i, j, k).y, u(i, j + 1, k).y);
            du_ds(i, j, k).z =
              slope_factor *
              L::limit(u(i, j - 1, k).z, u(i, j, k).z, u(i, j + 1, k).z);
            // Primitive quantities
            uTail(i, j, k) = u(i, j, k) + extrap_factor * du_ds(i, j, k);
            uHead(i, j, k) = u(i, j, k) - extrap_factor * du_ds(i, j, k);

            // Conserved quantities
            ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
            ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);

            // Compute slope
            dp_ds(i, j, k) =
              slope_factor *
              L::limit(p(i, j - 1, k), p(i, j, k), p(i, j + 1, k));

            // Primitive quantities
            pTail(i, j, k) = p(i, j, k) + extrap_factor * dp_ds(i, j, k);
            pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

            // Conserved quantities
            rETail(i, j, k) =
              p2rEFactor * pTail(i, j, k) +
              0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
            rEHead(i, j, k) =
              p2rEFactor * pHead(i, j, k) +
              0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();

            /*----------------------------------------------------------------------*
              Predictor step: approximate flux by averaging the governing flux
              functions at either side of the cell.
            *----------------------------------------------------------------------*/

            // Density
            q(i, j, k) =
              r(i, j, k) - courant * (ruTail(i, j, k).y - ruHead(i, j, k).y);

            // Momentum
            // ruv
            // rv^2+p
            // rvw
            const vec<3> f_qu_Tail{
              rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).y,
              ruTail(i, j, k).y * uTail(i, j, k).y + pTail(i, j, k),
              rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z};
            const vec<3> f_qu_Head{
              rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).y,
              ruHead(i, j, k).y * uHead(i, j, k).y + pHead(i, j, k),
              rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z};
            qu(i, j, k) = ru(i, j, k) - courant * (f_qu_Tail - f_qu_Head);

            // Total Energy
            const double f_qE_Tail =
              (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).y;
            const double f_qE_Head =
              (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).y;
            qE(i, j, k) = rE(i, j, k) - courant * (f_qE_Tail - f_qE_Head);
          } // for
        }
      };

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_primitives_3dy") {
        for(auto j : m.template cells<ax::y, dm::predictor>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // Primitives
            u(i, j, k) = qu(i, j, k) / q(i, j, k);
            p(i, j, k) =
              (gamma - 1.0) *
              (qE(i, j, k) - 0.5 * q(i, j, k) * u(i, j, k).norm_squared());
          } // for
        }
      };

        /*----------------------------------------------------------------------*
          Radiation slope, face, energy.
         *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(k,
        (m.template cells<ax::z, dm::quantities>()),
        "adv_slope_face_energy_rad_3dy") {
        for(auto j : m.template cells<ax::y, dm::predictor>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // Compute slope
            dErad_ds(i, j, k) =
              slope_factor *
              L::limit(Erad(i, j - 1, k), Erad(i, j, k), Erad(i, j + 1, k));

            // Conserved quantities
            EradTail(i, j, k) =
              Erad(i, j, k) + extrap_factor * dErad_ds(i, j, k);
            EradHead(i, j, k) =
              Erad(i, j, k) - extrap_factor * dErad_ds(i, j, k);

            // Radiation Energy
            qErad(i, j, k) =
              Erad(i, j, k) -
              courant * ((EradTail(i, j, k) * uTail(i, j, k).y) -
                          (EradHead(i, j, k) * uHead(i, j, k).y));
          } // for
        }
      };

#endif

      forall(k,
        (m.template cells<ax::z, dm::quantities>()),
        "adv_re_face_riemann_3dy") {
        for(auto j : m.template cells<ax::y, dm::corrector>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            /*----------------------------------------------------------------------*
              Reconstruct faces from intermediates.
            *----------------------------------------------------------------------*/
            // Primitive quantities
            rTail(i, j, k) =
              q(i, j - 1, k) + extrap_factor * dr_ds(i, j - 1, k);
            rHead(i, j, k) = q(i, j, k) - extrap_factor * dr_ds(i, j, k);
            uTail(i, j, k) =
              u(i, j - 1, k) + extrap_factor * du_ds(i, j - 1, k);
            uHead(i, j, k) = u(i, j, k) - extrap_factor * du_ds(i, j, k);
            pTail(i, j, k) =
              p(i, j - 1, k) + extrap_factor * dp_ds(i, j - 1, k);
            pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

            // Conserved quantities
            ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
            ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);
            rETail(i, j, k) =
              p2rEFactor * pTail(i, j, k) +
              0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
            rEHead(i, j, k) =
              p2rEFactor * pHead(i, j, k) +
              0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();
            EradTail(i, j, k) =
              qErad(i, j - 1, k) + extrap_factor * dErad_ds(i, j - 1, k);
            EradHead(i, j, k) =
              qErad(i, j, k) - extrap_factor * dErad_ds(i, j, k);

            /*----------------------------------------------------------------------*
              Compute corrector setp: Riemann solver.
            *----------------------------------------------------------------------*/

            // Update min/max eigenvalues on tail face.
            const double cT =
              std::sqrt(gamma * pTail(i, j, k) / rTail(i, j, k));
            const double LminT = uTail(i, j, k).y - cT;
            const double LmaxT = uTail(i, j, k).y + cT;

            // Update min/max eigenvalues on head face.
            const double cH =
              std::sqrt(gamma * pHead(i, j, k) / rHead(i, j, k));
            const double LminH = uHead(i, j, k).y - cH;
            const double LmaxH = uHead(i, j, k).y + cH;

            // Fluxes from left and right state, for the Riemann solver
            const double f_r_T{ruTail(i, j, k).y};
            const double f_r_H{ruHead(i, j, k).y};
            const vec<3> f_ru_T{
              rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).y,
              ruTail(i, j, k).y * uTail(i, j, k).y + pTail(i, j, k),
              rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z};
            const vec<3> f_ru_H{
              rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).y,
              ruHead(i, j, k).y * uHead(i, j, k).y + pHead(i, j, k),
              rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z};
            const double f_rE_T{
              (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).y};
            const double f_rE_H{
              (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).y};

            const auto hll_fluxes =
              numerical_algorithms::hll_hydro(rTail(i, j, k),
                rHead(i, j, k),
                ruTail(i, j, k),
                ruHead(i, j, k),
                rETail(i, j, k),
                rEHead(i, j, k),
                f_r_T,
                f_r_H,
                f_ru_T,
                f_ru_H,
                f_rE_T,
                f_rE_H,
                LminT,
                LmaxT,
                LminH,
                LmaxH);

            rF(i, j, k) = get<0>(hll_fluxes);
            ruF(i, j, k) = get<1>(hll_fluxes);
            rEF(i, j, k) = get<2>(hll_fluxes);
          } // for
        }
      };

      /*----------------------------------------------------------------------*
        Update averages (in-place).
       *----------------------------------------------------------------------*/

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_r_3dy") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            r(i, j, k) =
              r(i, j, k) + update_factor * (rF(i, j, k) - rF(i, j + 1, k));
          } // for
        }
      };

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_ru_3dy") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            ru(i, j, k) =
              ru(i, j, k) + update_factor * (ruF(i, j, k) - ruF(i, j + 1, k));
          } // for
        }
      };

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_rE_3dy") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            rE(i, j, k) =
              rE(i, j, k) + update_factor * (rEF(i, j, k) - rEF(i, j + 1, k));
          } // for
        }
      };

        /*----------------------------------------------------------------------*
          Radiation update (in-place).
         *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_re_face_rad_3dy") {
        for(auto j : m.template cells<ax::y, dm::corrector>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // Conserved quantities
            EradTail(i, j, k) =
              qErad(i, j - 1, k) + extrap_factor * dErad_ds(i, j - 1, k);
            EradHead(i, j, k) =
              qErad(i, j, k) - extrap_factor * dErad_ds(i, j, k);

            EradF(i, j, k) =
              numerical_algorithms::advect_Erad(EradTail(i, j, k),
                EradHead(i, j, k),
                uTail(i, j, k).y,
                uHead(i, j, k).y);
          } // for
        }
      };

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_Erad_3dy") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            Erad(i, j, k) =
              Erad(i, j, k) +
              update_factor * (EradF(i, j, k) - EradF(i, j + 1, k));
          } // for
        }
      };

#endif
    }
    else if constexpr(A == ax::z) {
      const double slope_factor{1.0 / m.template delta<A>()};
      const double extrap_factor{m.template delta<A>() / 2.0};
      const double p2rEFactor{1.0 / (gamma - 1.0)};
      const double courant{dt / (2.0 * m.template delta<A>())};

      forall(k,
        (m.template cells<ax::z, dm::predictor>()),
        "adv_slope_face_energy_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            /*----------------------------------------------------------------------*
              Compute slopes and reconstruct faces.
            *----------------------------------------------------------------------*/
            // Compute slope
            dr_ds(i, j, k) =
              slope_factor *
              L::limit(r(i, j, k - 1), r(i, j, k), r(i, j, k + 1));
            // Primitive quantities
            rTail(i, j, k) = r(i, j, k) + extrap_factor * dr_ds(i, j, k);
            rHead(i, j, k) = r(i, j, k) - extrap_factor * dr_ds(i, j, k);

            // Compute slope
            du_ds(i, j, k).x =
              slope_factor *
              L::limit(u(i, j, k - 1).x, u(i, j, k).x, u(i, j, k + 1).x);
            du_ds(i, j, k).y =
              slope_factor *
              L::limit(u(i, j, k - 1).y, u(i, j, k).y, u(i, j, k + 1).y);
            du_ds(i, j, k).z =
              slope_factor *
              L::limit(u(i, j, k - 1).z, u(i, j, k).z, u(i, j, k + 1).z);
            // Primitive quantities
            uTail(i, j, k) = u(i, j, k) + extrap_factor * du_ds(i, j, k);
            uHead(i, j, k) = u(i, j, k) - extrap_factor * du_ds(i, j, k);

            // Conserved quantities
            ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
            ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);

            // Compute slope
            dp_ds(i, j, k) =
              slope_factor *
              L::limit(p(i, j, k - 1), p(i, j, k), p(i, j, k + 1));
            // Primitive quantities
            pTail(i, j, k) = p(i, j, k) + extrap_factor * dp_ds(i, j, k);
            pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

            // Conserved quantities
            rETail(i, j, k) =
              p2rEFactor * pTail(i, j, k) +
              0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
            rEHead(i, j, k) =
              p2rEFactor * pHead(i, j, k) +
              0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();

            /*----------------------------------------------------------------------*
              Predictor step: approximate flux by averaging the governing flux
              functions at either side of the cell.
            *----------------------------------------------------------------------*/
            // Density
            q(i, j, k) =
              r(i, j, k) - courant * (ruTail(i, j, k).z - ruHead(i, j, k).z);

            // Momentum
            // ruw
            // rvw
            // rw^2+p
            const vec<3> f_qu_Tail{
              rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).z,
              rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z,
              ruTail(i, j, k).z * uTail(i, j, k).z + pTail(i, j, k)};
            const vec<3> f_qu_Head{
              rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).z,
              rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z,
              ruHead(i, j, k).z * uHead(i, j, k).z + pHead(i, j, k)};
            qu(i, j, k) = ru(i, j, k) - courant * (f_qu_Tail - f_qu_Head);

            // Total Energy
            const double f_qE_Tail =
              (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).z;
            const double f_qE_Head =
              (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).z;
            qE(i, j, k) = rE(i, j, k) - courant * (f_qE_Tail - f_qE_Head);
          } // for
        }
      };

      forall(
        k, (m.template cells<ax::z, dm::predictor>()), "adv_primitives_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // Primitives
            u(i, j, k) = qu(i, j, k) / q(i, j, k);
            p(i, j, k) =
              (gamma - 1.0) *
              (qE(i, j, k) - 0.5 * q(i, j, k) * u(i, j, k).norm_squared());
          } // for
        }
      };

        /*----------------------------------------------------------------------*
          Radiation slope, face, energy.
         *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(k,
        (m.template cells<ax::z, dm::predictor>()),
        "adv_slope_face_energy_rad_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // Slope
            dErad_ds(i, j, k) =
              slope_factor *
              L::limit(Erad(i, j, k - 1), Erad(i, j, k), Erad(i, j, k + 1));
            // Conserved quantities
            EradTail(i, j, k) =
              Erad(i, j, k) + extrap_factor * dErad_ds(i, j, k);
            EradHead(i, j, k) =
              Erad(i, j, k) - extrap_factor * dErad_ds(i, j, k);

            // Radiation energy
            qErad(i, j, k) =
              Erad(i, j, k) -
              courant * ((EradTail(i, j, k) * uTail(i, j, k).z) -
                          (EradHead(i, j, k) * uHead(i, j, k).z));
          } // for
        }
      };

#endif

      forall(k,
        (m.template cells<ax::z, dm::corrector>()),
        "adv_re_face_riemann_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            /*----------------------------------------------------------------------*
              Reconstruct faces from intermediates.
            *----------------------------------------------------------------------*/
            // Primitive quantities
            rTail(i, j, k) =
              q(i, j, k - 1) + extrap_factor * dr_ds(i, j, k - 1);
            rHead(i, j, k) = q(i, j, k) - extrap_factor * dr_ds(i, j, k);
            uTail(i, j, k) =
              u(i, j, k - 1) + extrap_factor * du_ds(i, j, k - 1);
            uHead(i, j, k) = u(i, j, k) - extrap_factor * du_ds(i, j, k);
            pTail(i, j, k) =
              p(i, j, k - 1) + extrap_factor * dp_ds(i, j, k - 1);
            pHead(i, j, k) = p(i, j, k) - extrap_factor * dp_ds(i, j, k);

            // Conserved quantities
            ruTail(i, j, k) = rTail(i, j, k) * uTail(i, j, k);
            ruHead(i, j, k) = rHead(i, j, k) * uHead(i, j, k);
            rETail(i, j, k) =
              p2rEFactor * pTail(i, j, k) +
              0.5 * rTail(i, j, k) * uTail(i, j, k).norm_squared();
            rEHead(i, j, k) =
              p2rEFactor * pHead(i, j, k) +
              0.5 * rHead(i, j, k) * uHead(i, j, k).norm_squared();

            /*----------------------------------------------------------------------*
              Compute corrector setp: Riemann solver.
            *----------------------------------------------------------------------*/
            // Update min/max eigenvalues on tail face.
            const double cT =
              std::sqrt(gamma * pTail(i, j, k) / rTail(i, j, k));
            const double LminT = uTail(i, j, k).z - cT;
            const double LmaxT = uTail(i, j, k).z + cT;

            // Update min/max eigenvalues on head face.
            const double cH =
              std::sqrt(gamma * pHead(i, j, k) / rHead(i, j, k));
            const double LminH = uHead(i, j, k).z - cH;
            const double LmaxH = uHead(i, j, k).z + cH;

            // Fluxes from left and right state, for the Riemann solver
            const double f_r_T{ruTail(i, j, k).z};
            const double f_r_H{ruHead(i, j, k).z};
            const vec<3> f_ru_T{
              rTail(i, j, k) * uTail(i, j, k).x * uTail(i, j, k).z,
              rTail(i, j, k) * uTail(i, j, k).y * uTail(i, j, k).z,
              ruTail(i, j, k).z * uTail(i, j, k).z + pTail(i, j, k)};
            const vec<3> f_ru_H{
              rHead(i, j, k) * uHead(i, j, k).x * uHead(i, j, k).z,
              rHead(i, j, k) * uHead(i, j, k).y * uHead(i, j, k).z,
              ruHead(i, j, k).z * uHead(i, j, k).z + pHead(i, j, k)};
            const double f_rE_T{
              (rETail(i, j, k) + pTail(i, j, k)) * uTail(i, j, k).z};
            const double f_rE_H{
              (rEHead(i, j, k) + pHead(i, j, k)) * uHead(i, j, k).z};

            const auto hll_fluxes =
              numerical_algorithms::hll_hydro(rTail(i, j, k),
                rHead(i, j, k),
                ruTail(i, j, k),
                ruHead(i, j, k),
                rETail(i, j, k),
                rEHead(i, j, k),
                f_r_T,
                f_r_H,
                f_ru_T,
                f_ru_H,
                f_rE_T,
                f_rE_H,
                LminT,
                LmaxT,
                LminH,
                LmaxH);

            rF(i, j, k) = get<0>(hll_fluxes);
            ruF(i, j, k) = get<1>(hll_fluxes);
            rEF(i, j, k) = get<2>(hll_fluxes);
          } // for
        }
      };

      /*----------------------------------------------------------------------*
        Update averages (in-place).
       *----------------------------------------------------------------------*/

      const double update_factor = dt / m.template delta<A>();

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_r_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            r(i, j, k) =
              r(i, j, k) + update_factor * (rF(i, j, k) - rF(i, j, k + 1));
          } // for
        }
      };

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_ru_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            ru(i, j, k) =
              ru(i, j, k) + update_factor * (ruF(i, j, k) - ruF(i, j, k + 1));
          } // for
        }
      };

      forall(k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_rE_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            rE(i, j, k) =
              rE(i, j, k) + update_factor * (rEF(i, j, k) - rEF(i, j, k + 1));
          } // for
        }
      };

        /*----------------------------------------------------------------------*
          Radiation update (in-place).
         *----------------------------------------------------------------------*/
#ifndef DISABLE_RADIATION

      forall(
        k, (m.template cells<ax::z, dm::corrector>()), "adv_advect_rad_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // Conserved quantities
            EradTail(i, j, k) =
              qErad(i, j, k - 1) + extrap_factor * dErad_ds(i, j, k - 1);
            EradHead(i, j, k) =
              qErad(i, j, k) - extrap_factor * dErad_ds(i, j, k);

            EradF(i, j, k) =
              numerical_algorithms::advect_Erad(EradTail(i, j, k),
                EradHead(i, j, k),
                uTail(i, j, k).z,
                uHead(i, j, k).z);
          } // for
        }
      };

      forall(
        k, (m.template cells<ax::z, dm::quantities>()), "adv_upd_Erad_3dz") {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            Erad(i, j, k) =
              Erad(i, j, k) +
              update_factor * (EradF(i, j, k) - EradF(i, j, k + 1));
          } // for
        }
      };

#endif
    } // if ax::z
  } // if ax==z
} // advance

} // namespace hard::tasks::hydro

#endif // HARD_TASKS_HYDRO_HH
