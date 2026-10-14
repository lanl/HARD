#ifndef HARD_HYDRO_MODULES_KELVIN_HELM_MM_HH
#define HARD_HYDRO_MODULES_KELVIN_HELM_MM_HH

namespace hard::tasks::initial_data {

struct kh_instability_mm_task : annotation_task_region {
  inline static const std::string name{"kh_instability_mm"};
};

/*----------------------------------------------------------------------------*
  Kelvin-Helmholtz Instability in MM.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
auto
kh_instability_mm(flecsi::exec::cpu,
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
  const std::vector<eos::eos_wrapper> & eos_mm) noexcept {

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

  const double wavenumber = 2.0 * M_PI;
  const double N = 4;

  // density, velocity and pressure for two fluids
  const double r_l = 1.0;
  const double u_l = -0.5;
  const double v_l = 0.0;
  const double p_l = 2.5;
  const int m_l = 0;

  const double r_h = 2.0;
  const double u_h = 0.5;
  const double v_h = 0.0;
  const double p_h = 2.5;
  const int m_h = 1;

  vec<D> vel;

  double sieL = 0., sieH = 0.;
  const double e_l =
    eos_mm[static_cast<std::size_t>(m_l)].ieRhoPSie(r_l, p_l, sieL);
  const double e_h =
    eos_mm[static_cast<std::size_t>(m_h)].ieRhoPSie(r_h, p_h, sieH);

  if constexpr(D == 1) {
    flog_fatal(
      "Kelvin-Helmholtz instability problem for D == 1 is not implemented")
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        // velocity perturbations in the Y-direction
        const double v_pert = 0.05 * sin(N * wavenumber * x);

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

          material(i, j) = m_l;
          frozen_material(i, j) = m_l;
          psi(i, j) = 0.25 - std::abs(y - 0.5);

          total_energy_density(i, j) =
            r_l * e_l + 0.5 * r_l * (vel.x() * vel.x() + vel.y() * vel.y());

          hard::tasks::util::pack_q<D>(q_scalar(i, j), r_l, u_l, p_l, e_l);
          q_scalar_known(i, j) = q_scalar(i, j);
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

          material(i, j) = m_h;
          frozen_material(i, j) = m_h;
          psi(i, j) = 0.25 - std::abs(y - 0.5);

          total_energy_density(i, j) =
            r_h * e_h + 0.5 * r_h * (vel.x() * vel.x() + vel.y() * vel.y());

          hard::tasks::util::pack_q<D>(q_scalar(i, j), r_h, u_h, p_h, e_h);
          q_scalar_known(i, j) = q_scalar(i, j);
        } // if
        // Initial exact constant extrapolation fields.
        // plus: left/minus material extended into plus/right side.
        hard::tasks::util::pack_q<D>(q_ext_plus(i, j), r_l, vel, p_l, e_l);

        // minus: right/plus material extended into minus/left side.
        hard::tasks::util::pack_q<D>(q_ext_minus(i, j), r_h, vel, p_h, e_h);
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
