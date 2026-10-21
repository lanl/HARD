#ifndef HARD_HYDRO_MODULES_KELVIN_HELM_HH
#define HARD_HYDRO_MODULES_KELVIN_HELM_HH

namespace hard::tasks::initial_data {

struct kh_instability_task : annotation_task_region {
  inline static const std::string name{"kh_instability"};
};

/*----------------------------------------------------------------------------*
  Kelvin-Helmholtz Instability.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
kh_instability(flecsi::exec::cpu,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  const eos::eos_wrapper & eos) noexcept {
  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  const double wavenumber = 2.0 * M_PI;
  const double N = 4;

  // density, velocity and pressure for two fluids
  const double r_l = 1.0;
  const double u_l = -0.5;
  const double v_l = 0.0;
  const double pL = 2.5;

  const double r_h = 2.0;
  const double u_h = 0.5;
  const double v_h = 0.0;
  const double pH = 2.5;

  vec<D> vel;

  double sieL = 0., sieH = 0.;
  const double e_l = eos.ieRhoPSie(r_l, pL, sieL);
  const double e_h = eos.ieRhoPSie(r_h, pH, sieH);

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        double v_pert;
        // velocity perturbations in the Y-direction
        v_pert = 0.05 * sin(N * wavenumber * x);

        // initialize two different density and velocity fluids
        if(std::abs(y - 0.5) > 0.25) {
          mass_density(i, j) = r_l;
          vel.x() = u_l;
          momentum_density(i, j).x() = r_l * vel.x();
          if(std::abs(y - 0.25) < 0.1 || std::abs(y - 0.75) < 0.1) {
            vel.y() = v_pert;
          }
          else {
            vel.y() = v_l;
          }
          momentum_density(i, j).y() = r_l * vel.y();

          total_energy_density(i, j) =
            r_l * e_l + 0.5 * r_l * (vel.x() * vel.x() + vel.y() * vel.y());
        }
        else {
          mass_density(i, j) = r_h;
          vel.x() = u_h;
          momentum_density(i, j).x() = r_h * vel.x();
          if(std::abs(y - 0.25) < 0.1 || std::abs(y - 0.75) < 0.1) {
            vel.y() = v_pert;
          }
          else {
            vel.y() = v_h;
          }
          momentum_density(i, j).y() = r_h * vel.y();

          total_energy_density(i, j) =
            r_h * e_h + 0.5 * r_h * (vel.x() * vel.x() + vel.y() * vel.y());
        } // if
      } // for
    }; // forall
  }
  else /* D == 3 */ {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 3 is not implemented")
  } // if
} //  kh_instability

} // namespace hard::tasks::initial_data

#endif // HARD_HYDRO_MODULES_KELVIN_HELM_HH
