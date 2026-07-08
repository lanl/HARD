#ifndef HARD_HYDRO_MODULES_KELVIN_HELM_MM_HH
#define HARD_HYDRO_MODULES_KELVIN_HELM_MM_HH

namespace hard::tasks::initial_data {

/*----------------------------------------------------------------------------*
  Kelvin-Helmholtz Instability.
 *----------------------------------------------------------------------------*/
// !!! NOT INITIALIZED OR LINKED YET !!!

template<std::size_t D>
auto
kh_instability_mm(flecsi::exec::cpu s,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<D>>::template accessor<wo, na> momentum_density_a,
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
  const std::vector<eos::eos_wrapper> & eos_mm) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
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

  // setting parameters of ghost cells
  auto pack_q = [](auto & qcell,
                  double rho,
                  const vec<D> & velocity,
                  double pressure,
                  double sie) {
    qcell.rho() = rho;
    qcell.u() = velocity;
    qcell.p() = pressure;
    qcell.e() = sie;
    qcell.n_comps() = 0;
  };

  const double wavenumber = 2.0 * M_PI;
  const double N = 4;

  // density, velocity and pressure for two fluids
  const double rL = 1.0;
  const double uL = -0.5;
  const double vL = 0.0;
  const double pL = 2.5;
  const int mL = 0;

  const double rH = 2.0;
  const double uH = 0.5;
  const double vH = 0.0;
  const double pH = 2.5;
  const int mH = 1;

  // vec<D> uL, uH; // unsure
  const double eL = util::find_sie(eos_mm[mL], rL, pL);
  const double eH = util::find_sie(eos_mm[mH], rH, pH);

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        // initialize two different density and velocity fluids
        if(std::abs(y - 0.5) > 0.25) {
          mass_density(i, j) = rL;
          momentum_density(i, j).x() = rL * uL;
          momentum_density(i, j).y() = rL * vL;

          material(i, j) = mL;
          frozen_material(i, j) = mL;
          psi(i, j) = std::abs(y - 0.5); // is this right ?

          const double e =
            util::find_sie(eos_mm[static_cast<std::size_t>(mL)], rL, pL);
          total_energy_density(i, j) = rL * e + 0.5 * rL * (vL * vL);

          pack_q(q_scalar(i, j), rL, uL, pL, eL);
          q_scalar_known(i, j) = q_scalar(i, j);
        }
        else {
          mass_density(i, j) = rH;
          momentum_density(i, j).x() = rH * uH;
          momentum_density(i, j).y() = rH * vH;

          material(i, j) = mH;
          frozen_material(i, j) = mH;
          psi(i, j) = std::abs(y - 0.5); // is this right ?

          const double e =
            util::find_sie(eos_mm[static_cast<std::size_t>(mH)], rH, pH);
          total_energy_density(i, j) = rH * e + 0.5 * rH * (vH * vH);

          pack_q(q_scalar(i, j), rH, uH, pH, eH);
          q_scalar_known(i, j) = q_scalar(i, j);
        } // if

        // velocity perturbations in the Y-direction
        if(std::abs(y - 0.25) < 0.1) {
          momentum_density(i, j).y() = 0.05 * sin(N * wavenumber * x);
        }
        if(std::abs(y - 0.75) < 0.1) {
          momentum_density(i, j).y() = 0.05 * sin(N * wavenumber * x);
        }
      } // for
    }; // forall
  }
  else /* D == 3 */ {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 3 is not implemented")
  } // if
} //  kh_instability

} // namespace hard::tasks::initial_data

#endif // HARD_HYDRO_MODULES_KELVIN_HELM_MM_HH