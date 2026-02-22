#ifndef HARD_MODULES_HYDRO_TASKS_OPERATORS_HH
#define HARD_MODULES_HYDRO_TASKS_OPERATORS_HH

namespace hard::tasks::rad {

  
template<std::size_t D>
void
full_weighting(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> mf,
  typename mesh<D>::template accessor<ro> mc,
  typename field<double>::template accessor<ro, ro> rfa,
  typename field<double>::template accessor<wo, ro> fca) noexcept {
  // TODO: fca could be <wo, na>, since only writing quantities

  auto rf = mf.template mdcolex<is::cells>(rfa);
  auto fc = mc.template mdcolex<is::cells>(fca);

  if constexpr(D == 1) {
    s.executor().forall(i, (mc.template cells<ax::x, dm::quantities>())) {
      auto fi = 2 * i;
      fc(i) = 0.125 * (3.0 * (rf(fi - 2) + rf(fi - 1)) + rf(fi - 3) + rf(fi));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(fc,
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      auto fj = 2 * j;
      auto fi = 2 * i;
      fc(i, j) =
        0.015625 *
        (9.0 * (rf(fi - 2, fj - 2) + rf(fi - 1, fj - 2) + rf(fi - 2, fj - 1) +
                 rf(fi - 1, fj - 1)) +

          3.0 * (rf(fi - 3, fj - 2) + rf(fi - 3, fj - 1) + rf(fi, fj - 2) +
                  rf(fi, fj - 1) + rf(fi - 2, fj - 3) + rf(fi - 1, fj - 3) +
                  rf(fi - 2, fj) + rf(fi - 1, fj)) +

          rf(fi, fj - 3) + rf(fi - 3, fj - 3) + rf(fi - 3, fj) + rf(fi, fj));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(fc,
      mc.template cells<ax::z, dm::quantities>(),
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      auto fk = 2 * k;
      auto fj = 2 * j;
      auto fi = 2 * i;

      fc(i, j, k) =
        0.001953125 *
        (27.0 * (rf(fi - 2, fj - 2, fk - 2) + rf(fi - 2, fj - 2, fk - 1) +
                  rf(fi - 2, fj - 1, fk - 2) + rf(fi - 1, fj - 2, fk - 2) +
                  rf(fi - 2, fj - 1, fk - 1) + rf(fi - 1, fj - 2, fk - 1) +
                  rf(fi - 1, fj - 1, fk - 2) + rf(fi - 1, fj - 1, fk - 1)) +

          9.0 * (rf(fi - 2, fj - 2, fk - 3) + rf(fi - 2, fj - 3, fk - 2) +
                  rf(fi - 3, fj - 2, fk - 2) + rf(fi - 2, fj - 1, fk - 3) +
                  rf(fi - 2, fj - 3, fk - 1) + rf(fi - 3, fj - 2, fk - 1) +
                  rf(fi - 1, fj - 2, fk - 3) + rf(fi - 1, fj - 3, fk - 2) +
                  rf(fi - 3, fj - 1, fk - 2) + rf(fi - 1, fj - 1, fk - 3) +
                  rf(fi - 1, fj - 3, fk - 1) + rf(fi - 3, fj - 1, fk - 1) +

                  rf(fi - 2, fj - 2, fk) + rf(fi - 2, fj, fk - 2) +
                  rf(fi, fj - 2, fk - 2) + rf(fi - 2, fj - 1, fk) +
                  rf(fi - 2, fj, fk - 1) + rf(fi, fj - 2, fk - 1) +
                  rf(fi - 1, fj - 2, fk) + rf(fi - 1, fj, fk - 2) +
                  rf(fi, fj - 1, fk - 2) + rf(fi - 1, fj - 1, fk) +
                  rf(fi - 1, fj, fk - 1) + rf(fi, fj - 1, fk - 1)) +

          3.0 *
            (rf(fi - 2, fj - 3, fk - 3) + rf(fi - 3, fj - 2, fk - 3) +
              rf(fi - 3, fj - 3, fk - 2) + rf(fi - 1, fj - 3, fk - 3) +
              rf(fi - 3, fj - 1, fk - 3) + rf(fi - 3, fj - 3, fk - 1) +

              rf(fi - 2, fj, fk - 3) + rf(fi, fj - 2, fk - 3) +
              rf(fi, fj - 3, fk - 2) + rf(fi - 1, fj, fk - 3) +
              rf(fi, fj - 1, fk - 3) + rf(fi, fj - 3, fk - 1) +

              rf(fi - 2, fj - 3, fk) + rf(fi - 3, fj - 2, fk) +
              rf(fi - 3, fj, fk - 2) + rf(fi - 1, fj - 3, fk) +
              rf(fi - 3, fj - 1, fk) + rf(fi - 3, fj, fk - 1) +

              rf(fi - 2, fj, fk) + rf(fi, fj - 2, fk) + rf(fi, fj, fk - 2) +
              rf(fi - 1, fj, fk) + rf(fi, fj - 1, fk) + rf(fi, fj, fk - 1)) +

          rf(fi - 3, fj - 3, fk - 3) + rf(fi - 3, fj - 3, fk) +
          rf(fi - 3, fj, fk - 3) + rf(fi, fj - 3, fk - 3) + rf(fi - 3, fj, fk) +
          rf(fi, fj - 3, fk) + rf(fi, fj, fk - 3) + rf(fi, fj, fk));
    }; // for
  } // if
} // full_weighting

template<std::size_t D>
void
nlinear_interpolation(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> mc,
  typename mesh<D>::template accessor<ro> mf,
  typename field<double>::template accessor<ro, ro> cfa,
  typename field<double>::template accessor<wo, ro> ffa) noexcept {
  // TODO: As above, ffa could be <wo, na>, since only writing quantities

  auto cf = mc.template mdcolex<is::cells>(cfa);
  auto ff = mf.template mdcolex<is::cells>(ffa);

  if constexpr(D == 1) {
    s.executor().forall(i, (mf.template cells<ax::x, dm::quantities>())) {
      auto ci = i / 2 + 1;
      auto di = +1;
      if((i + 1) % 2) {
        di = -1;
      }
      ff(i) = 0.25 * (3 * cf(ci) + cf(ci + di));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(cf,
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      auto cj = j / 2 + 1;
      auto dj = 1;
      if((j + 1) % 2) {
        dj = -1;
      }
      auto ci = i / 2 + 1;
      auto di = 1;
      if((i + 1) % 2) {
        di = -1;
      }
      ff(i, j) = 0.0625 * (9 * cf(ci, cj) + 3 * cf(ci + di, cj) +
                            3 * cf(ci, cj + dj) + cf(ci + di, cj + dj));
    }; // for
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(cf,
      mc.template cells<ax::z, dm::quantities>(),
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      auto ck = k / 2 + 1;
      auto dk = 1;
      if((k + 1) % 2) {
        dk = -1;
      }
      auto cj = j / 2 + 1;
      auto dj = 1;
      if((j + 1) % 2) {
        dj = -1;
      }
      auto ci = i / 2 + 1;
      auto di = 1;
      if((i + 1) % 2) {
        di = -1;
      }
      ff(i, j, k) =
        0.015625 *
        (27 * cf(ci, cj, ck) + 9 * cf(ci + di, cj, ck) +
          9 * cf(ci, cj + dj, ck) + 9 * cf(ci, cj, ck + dk) +
          3 * cf(ci + di, cj + dj, ck) + 3 * cf(ci + di, cj, ck + dk) +
          3 * cf(ci, cj + dj, ck + dk) + cf(ci + di, cj + dj, ck + dk));
    };
  } // if
} // nlinear_interpolation

template<std::size_t D>
void
cell_centered_weighting(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> mf,
  typename mesh<D>::template accessor<ro> mc,
  field<double>::accessor<ro, na> rfa,
  field<double>::accessor<wo, na> fca) noexcept {

  auto rf = mf.template mdcolex<is::cells>(rfa);
  auto fc = mc.template mdcolex<is::cells>(fca);

  if constexpr(D == 1) {
    s.executor().forall(i, (mc.template cells<ax::x, dm::quantities>())) {
      auto fi = 2 * i - mc.ghost_zone_size();
      fc(i) = 0.5 * (rf(fi) + rf(fi + 1));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(fc,
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      auto fj = 2 * j - mc.ghost_zone_size();
      auto fi = 2 * i - mc.ghost_zone_size();

      fc(i, j) = 0.25 * (rf(fi, fj) + rf(fi + 1, fj) + rf(fi, fj + 1) +
                          rf(fi + 1, fj + 1));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(fc,
      mc.template cells<ax::z, dm::quantities>(),
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      auto fk = 2 * k - mc.ghost_zone_size();
      auto fj = 2 * j - mc.ghost_zone_size();
      auto fi = 2 * i - mc.ghost_zone_size();

      fc(i, j, k) =
        0.125 *
        (rf(fi, fj, fk) + rf(fi + 1, fj, fk) + rf(fi, fj + 1, fk) +
          rf(fi, fj, fk + 1) + rf(fi + 1, fj + 1, fk) + rf(fi + 1, fj, fk + 1) +
          rf(fi, fj + 1, fk + 1) + rf(fi + 1, fj + 1, fk + 1));

    }; // for
  } // if
} // full_weighting

template<std::size_t D>
void
cell_centered_interpolation(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> mc,
  typename mesh<D>::template accessor<ro> mf,
  field<double>::accessor<ro, na> cfa,
  field<double>::accessor<wo, na> ffa) noexcept {

  auto cf = mc.template mdcolex<is::cells>(cfa);
  auto ff = mf.template mdcolex<is::cells>(ffa);

  if constexpr(D == 1) {
    s.executor().forall(i, (mc.template cells<ax::x, dm::quantities>())) {
      auto fi = 2 * i - mc.ghost_zone_size();
      ff(i) = cf(fi);
      ff(i + 1) = cf(fi);
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(cf,
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      auto fj = 2 * j - mc.ghost_zone_size();
      auto fi = 2 * i - mc.ghost_zone_size();

      ff(fi, fj) = cf(i, j);
      ff(fi + 1, fj) = cf(i, j);
      ff(fi, fj + 1) = cf(i, j);
      ff(fi + 1, fj + 1) = cf(i, j);

    }; // for
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(cf,
      mc.template cells<ax::z, dm::quantities>(),
      mc.template cells<ax::y, dm::quantities>(),
      mc.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      auto fj = 2 * j - mc.ghost_zone_size();
      auto fi = 2 * i - mc.ghost_zone_size();
      auto fk = 2 * k - mc.ghost_zone_size();

      ff(fi, fj, fk) = cf(i, j, k);
      ff(fi + 1, fj, fk) = cf(i, j, k);
      ff(fi, fj + 1, fk) = cf(i, j, k);
      ff(fi, fj, fk + 1) = cf(i, j, k);
      ff(fi + 1, fj + 1, fk) = cf(i, j, k);
      ff(fi + 1, fj + 1, fk + 1) = cf(i, j, k);
      ff(fi + 1, fj, fk + 1) = cf(i, j, k);
      ff(fi, fj + 1, fk + 1) = cf(i, j, k);
    };
  } // if
}
  
template<std::size_t D>
void
apply_operator(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  field<double>::accessor<wo, ro> ua_new,
  field<double>::accessor<ro, ro> ua_old) noexcept {

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u_new = m.template mdcolex<is::cells>(ua_new);
  auto u_old = m.template mdcolex<is::cells>(ua_old);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      u_new(i) = (Ew(i)[dirs::c] * u_old(i) - Ew(i)[dirs::w] * u_old(i - 1) -
                  Ew(i + 1)[dirs::w] * u_old(i + 1));
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u_new,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      u_new(i, j) =
        (Ew(i, j)[dirs::c] * u_old(i, j) - Ew(i, j)[dirs::w] * u_old(i - 1, j) -
          Ew(i + 1, j)[dirs::w] * u_old(i + 1, j) -
          Ew(i, j)[dirs::s] * u_old(i, j - 1) -
          Ew(i, j + 1)[dirs::s] * u_old(i, j + 1));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u_new,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      u_new(i, j, k) = (Ew(i, j, k)[dirs::c] * u_old(i, j, k) -
                        Ew(i, j, k)[dirs::w] * u_old(i - 1, j, k) -
                        Ew(i + 1, j, k)[dirs::w] * u_old(i + 1, j, k) -
                        Ew(i, j, k)[dirs::s] * u_old(i, j - 1, k) -
                        Ew(i, j + 1, k)[dirs::s] * u_old(i, j + 1, k) -
                        Ew(i, j, k)[dirs::d] * u_old(i, j, k - 1) -
                        Ew(i, j, k + 1)[dirs::d] * u_old(i, j, k + 1));
    }; // forall
  } // if
} // Ax_op


template<std::size_t D>
void
damped_jacobi(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<wo, na> ua_new,
  typename field<double>::template accessor<ro, ro> ua_old,
  typename field<double>::template accessor<ro, ro> fa,
  double omega) noexcept {

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u_new = m.template mdcolex<is::cells>(ua_new);
  auto u_old = m.template mdcolex<is::cells>(ua_old);
  auto f = m.template mdcolex<is::cells>(fa);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      const double z = (Ew(i)[dirs::w] * u_old(i - 1) +
                         Ew(i + 1)[dirs::w] * u_old(i + 1) + f(i)) /
                       Ew(i)[dirs::c];

      u_new(i) = u_old(i) + omega * (z - u_old(i));
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u_new,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      const double z = (Ew(i, j)[dirs::w] * u_old(i - 1, j) +
                         Ew(i + 1, j)[dirs::w] * u_old(i + 1, j) +
                         Ew(i, j)[dirs::s] * u_old(i, j - 1) +
                         Ew(i, j + 1)[dirs::s] * u_old(i, j + 1) + f(i, j)) /
                       Ew(i, j)[dirs::c];

      u_new(i, j) = u_old(i, j) + omega * (z - u_old(i, j));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u_new,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      const double z =
        (Ew(i, j, k)[dirs::w] * u_old(i - 1, j, k) +
          Ew(i + 1, j, k)[dirs::w] * u_old(i + 1, j, k) +
          Ew(i, j, k)[dirs::s] * u_old(i, j - 1, k) +
          Ew(i, j + 1, k)[dirs::s] * u_old(i, j + 1, k) +
          Ew(i, j, k)[dirs::d] * u_old(i, j, k - 1) +
          Ew(i, j, k + 1)[dirs::d] * u_old(i, j, k + 1) + f(i, j, k)) /
        Ew(i, j, k)[dirs::c];

      u_new(i, j, k) = u_old(i, j, k) + omega * (z - u_old(i, j, k));
    }; // forall
  } // if
} // damped_jacobi

  
template<std::size_t D>
void
residual(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<stencil<D>>::template accessor<ro, ro> Ew_a,
  typename field<double>::template accessor<ro, ro> ua,
  typename field<double>::template accessor<ro, ro> fa,
  typename field<double>::template accessor<wo, ro> ra) noexcept {
  // TODO: Not using the ghost cells of fa here, are we?

  auto Ew = m.template mdcolex<is::cells>(Ew_a);
  auto u = m.template mdcolex<is::cells>(ua);
  auto f = m.template mdcolex<is::cells>(fa);
  auto r = m.template mdcolex<is::cells>(ra);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      r(i) = f(i) - (Ew(i)[dirs::c] * u(i) - Ew(i)[dirs::w] * u(i - 1) -
                      Ew(i + 1)[dirs::w] * u(i + 1));
    };
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      r(i, j) = f(i, j) -
                (Ew(i, j)[dirs::c] * u(i, j) - Ew(i, j)[dirs::w] * u(i - 1, j) -
                  Ew(i + 1, j)[dirs::w] * u(i + 1, j) -
                  Ew(i, j)[dirs::s] * u(i, j - 1) -
                  Ew(i, j + 1)[dirs::s] * u(i, j + 1));
    }; // forall
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      r(i, j, k) = f(i, j, k) - (Ew(i, j, k)[dirs::c] * u(i, j, k) -
                                  Ew(i, j, k)[dirs::w] * u(i - 1, j, k) -
                                  Ew(i + 1, j, k)[dirs::w] * u(i + 1, j, k) -
                                  Ew(i, j, k)[dirs::s] * u(i, j - 1, k) -
                                  Ew(i, j + 1, k)[dirs::s] * u(i, j + 1, k) -
                                  Ew(i, j, k)[dirs::d] * u(i, j, k - 1) -
                                  Ew(i, j, k + 1)[dirs::d] * u(i, j, k + 1));
    }; // forall
  } // if
} // residual

template<std::size_t D>
void
correction(flecsi::exec::accelerator s,
  typename mesh<D>::template accessor<ro> m,
  typename field<double>::template accessor<rw, ro> ua,
  typename field<double>::template accessor<wo, ro> ea) noexcept {
  // TODO: Looks like ea should be <ro, na>

  auto u = m.template mdcolex<is::cells>(ua);
  auto e = m.template mdcolex<is::cells>(ea);

  if constexpr(D == 1) {
    s.executor().forall(i, (m.template cells<ax::x, dm::quantities>())) {
      u(i) += e(i);
    }; // for
  }
  else if constexpr(D == 2) {
    auto mdpolicy_qq = get_mdiota_policy(u,
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(ji, mdpolicy_qq) {
      auto [j, i] = ji;
      u(i, j) += e(i, j);
    }; // for
  }
  else /* D == 3 */ {
    auto mdpolicy_qqq = get_mdiota_policy(u,
      m.template cells<ax::z, dm::quantities>(),
      m.template cells<ax::y, dm::quantities>(),
      m.template cells<ax::x, dm::quantities>());

    s.executor().forall(kji, mdpolicy_qqq) {
      auto [k, j, i] = kji;
      u(i, j, k) += e(i, j, k);
    }; // for
  } // if
} // correction

} // namespace hard::tasks::rad

#endif // HARD_MODULE_RAD_TASKS_RAD_HH
