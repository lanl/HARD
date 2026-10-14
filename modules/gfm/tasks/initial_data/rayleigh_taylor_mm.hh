#ifndef HARD_HYDRO_MODULES_RAYLEIGH_TAYLOR_MM_HH
#define HARD_HYDRO_MODULES_RAYLEIGH_TAYLOR_MM_HH

namespace hard::tasks::initial_data {

struct rt_instability_mm_task : annotation_task_region {
  inline static const std::string name{"rt_instability_mm"};
};

/*----------------------------------------------------------------------------*
  Rayleigh-Taylor Instability in MM.
  Reference: https://www.astro.princeton.edu/~jstone/Athena/tests/rt/rt.html
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
rt_instability_mm(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  typename field<vec<D>>::template accessor<wo, na> gravity_accel_a,
  typename single<vec<D>>::template accessor<ro> gravity_const,
  field<double>::accessor<wo, na> total_energy_density_a,
  typename field<int>::accessor<wo, na> material_a,
  typename field<int>::accessor<wo, na> frozen_material_a,
  field<double>::accessor<wo, na> psi_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_scalar_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_scalar_known_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_ext_plus_a,
  typename field<hard::tasks::util::extp_q<D>>::template accessor<wo, na>
    q_ext_minus_a,
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto ga = m.template mdcolex<is::cells>(gravity_accel_a);
  auto g = gravity_const.get();
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto material = m.template mdcolex<is::cells>(material_a);
  auto frozen_material = m.template mdcolex<is::cells>(frozen_material_a);
  auto psi = m.template mdcolex<is::cells>(psi_a);
  auto q_scalar = m.template mdcolex<is::cells>(q_scalar_a);
  auto q_scalar_known = m.template mdcolex<is::cells>(q_scalar_known_a);
  auto q_ext_plus = m.template mdcolex<is::cells>(q_ext_plus_a);
  auto q_ext_minus = m.template mdcolex<is::cells>(q_ext_minus_a);

  // initial conditions
  const double p0 = 2.5;
  // For hydrostatic equilibrium pressure, add g*rho*y to p0
  const double v0 = 0.1;
  // For perturbed velocity, multiply by
  // ((1 + cos (4*pi*x))/2)*((1 + cos (3*pi*y))/2)

  // H is the upper fluid
  const double rH = 2.0;
  const double uH = 0.0;
  const int mH = 1;

  // L is the bottom fluid
  const double rL = 1.0;
  const double uL = 0.0;
  const int mL = 0;

  vec<D> vel;

  double sieL = 0., sieH = 0.;

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const auto x = m.template center<ax::x>(i);

      constexpr double x0 = 0.75; // interface of fluids
      const double v = v0 * 1 / 2 * std::exp(-std::pow((x - x0) / 0.05, 2));

      if(x >= x0) {
        // upper fluid
        mass_density(i) = rH;
        vel.x() = uH;
        ga(i).x() = rH * g.x();
        momentum_density(i).x() = rH * vel.x();
        const double pH = p0 + rL * g.x() * x0 + g.x() * rH * (x - x0);

        material(i) = mH;
        frozen_material(i) = mH;
        psi(i) = -std::abs(x - x0);

        const double eH =
          eos_mm[static_cast<std::size_t>(mH)].ieRhoPSie(rH, pH, sieH);
        total_energy_density(i) =
          rH * eH + 0.5 * rH * (vel.x() * vel.x() + v * v);

        hard::tasks::util::pack_q(q_scalar(i), rH, vel, pH, eH);
        q_scalar_known(i) = q_scalar(i);
        // Extrapolated states across the interface
        // plus: left/minus material extended into plus/right side.
        hard::tasks::util::pack_q(q_ext_plus(i), rH, vel, pH, eH);
      }
      else {
        // lower fluid
        mass_density(i) = rL;
        vel.x() = uL;
        ga(i).x() = g.x();
        momentum_density(i).x() = rL * vel.x();
        const double pL = p0 + rL * g.x() * x;

        material(i) = mL;
        frozen_material(i) = mL;
        psi(i) = -std::abs(x - x0);

        const double eL =
          eos_mm[static_cast<std::size_t>(mL)].ieRhoPSie(rL, pL, sieL);
        total_energy_density(i) =
          rL * eL + 0.5 * rL * (vel.x() * vel.x() + v * v);

        hard::tasks::util::pack_q(q_scalar(i), rL, vel, pL, eL);
        q_scalar_known(i) = q_scalar(i);
        // Extrapolated states across the interface
        // minus: right/plus material extended into minus/left side.
        hard::tasks::util::pack_q(q_ext_minus(i), rL, vel, pL, eL);
      } // if
    } // for
  } // D=1
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        constexpr double x0 = 0.25; // perturbation x location
        constexpr double y0 = 0.75; // interface of fluids
        const double v = v0 * (1 + cos(4 * M_PI * (x - x0))) / 2 *
                         std::exp(-std::pow((y - y0) / 0.05, 2));

        if(y >= y0) {
          // upper fluid
          mass_density(i, j) = rH;
          vel.x() = uH;
          vel.y() = v;
          ga(i, j).x() = g.x();
          ga(i, j).y() = g.y();
          momentum_density(i, j).x() = rH * vel.x();
          momentum_density(i, j).y() = rH * vel.y();
          const double pH = p0 + rL * g.y() * y0 + g.y() * rH * (y - y0);

          material(i, j) = mH;
          frozen_material(i, j) = mH;
          psi(i, j) = -std::abs(y - y0);

          const double eH =
            eos_mm[static_cast<std::size_t>(mH)].ieRhoPSie(rH, pH, sieH);
          total_energy_density(i, j) =
            rH * eH + 0.5 * rH * (vel.x() * vel.x() + vel.y() * vel.y());

          hard::tasks::util::pack_q(q_scalar(i, j), rH, vel, pH, eH);
          q_scalar_known(i, j) = q_scalar(i, j);
          // Extrapolated states across the interface
          // plus: left/minus material extended into plus/right side.
          hard::tasks::util::pack_q(q_ext_plus(i, j), rH, vel, pH, eH);
        }
        else {
          // lower fluid
          mass_density(i, j) = rL;
          vel.x() = uL;
          vel.y() = v;
          ga(i, j).x() = g.x();
          ga(i, j).y() = g.y();
          momentum_density(i, j).x() = rL * vel.x();
          momentum_density(i, j).y() = rL * vel.y();
          const double pL = p0 + rL * g.y() * y;

          material(i, j) = mL;
          frozen_material(i, j) = mL;
          psi(i, j) = -std::abs(y - y0);

          const double eL =
            eos_mm[static_cast<std::size_t>(mL)].ieRhoPSie(rL, pL, sieL);
          total_energy_density(i, j) =
            rL * eL + 0.5 * rL * (vel.x() * vel.x() + vel.y() * vel.y());

          hard::tasks::util::pack_q(q_scalar(i, j), rL, vel, pL, eL);
          q_scalar_known(i, j) = q_scalar(i, j);
          // Extrapolated states across the interface
          // minus: right/plus material extended into minus/left side.
          hard::tasks::util::pack_q(q_ext_minus(i, j), rL, vel, pL, eL);
        } // if
      } // for i
    } // for j
  } // D=2
  else {
    flog_fatal(
      "Rayleigh-Taylor instability problem for D == 3 is not implemented")
  } // D=3

} //  rt_instability_mm

} // namespace hard::tasks::initial_data

#endif // HARD_HYDRO_MODULES_RAYLEIGH_TAYLOR_MM_HH
