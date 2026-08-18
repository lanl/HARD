#ifndef HARD_MODULES_HYDRO_FOUR_SHOCK_JET_HH
#define HARD_MODULES_HYDRO_FOUR_SHOCK_JET_HH

#include "spec/types.hh"

namespace hard::tasks::initial_data {

//
// A four shock jet set-up
//
template<std::size_t Dim>
auto
four_shock_jet(flecsi::exec::cpu s,
  typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  const eos::eos_wrapper & eos) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);

  spec::config_py config(opt::config.value());

  const double xc{config["problem_parameters"]["xc"].cast<double>()};
  const double yc{config["problem_parameters"]["yc"].cast<double>()};

  struct values {
    double bot_left, top_left, top_right, bot_right;
  };

  const values rho_vals{
    .bot_left =
      config["problem_parameters"]["regions"]["bot_left"]["r"].cast<double>(),
    .top_left =
      config["problem_parameters"]["regions"]["top_left"]["r"].cast<double>(),
    .top_right =
      config["problem_parameters"]["regions"]["top_right"]["r"].cast<double>(),
    .bot_right =
      config["problem_parameters"]["regions"]["bot_right"]["r"].cast<double>(),
  };

  const values p_vals{
    .bot_left =
      config["problem_parameters"]["regions"]["bot_left"]["p"].cast<double>(),
    .top_left =
      config["problem_parameters"]["regions"]["top_left"]["p"].cast<double>(),
    .top_right =
      config["problem_parameters"]["regions"]["top_right"]["p"].cast<double>(),
    .bot_right =
      config["problem_parameters"]["regions"]["bot_right"]["p"].cast<double>(),
  };

  const values ux_vals{
    .bot_left =
      config["problem_parameters"]["regions"]["bot_left"]["ux"].cast<double>(),
    .top_left =
      config["problem_parameters"]["regions"]["top_left"]["ux"].cast<double>(),
    .top_right =
      config["problem_parameters"]["regions"]["top_right"]["ux"].cast<double>(),
    .bot_right =
      config["problem_parameters"]["regions"]["bot_right"]["ux"].cast<double>(),
  };

  const values uy_vals{
    .bot_left =
      config["problem_parameters"]["regions"]["bot_left"]["uy"].cast<double>(),
    .top_left =
      config["problem_parameters"]["regions"]["top_left"]["uy"].cast<double>(),
    .top_right =
      config["problem_parameters"]["regions"]["top_right"]["uy"].cast<double>(),
    .bot_right =
      config["problem_parameters"]["regions"]["bot_right"]["uy"].cast<double>(),
  };

  struct limits {
    double x0, x1, y0, y1;
  };

  const auto find_vals = [xc, yc](const values vals, const limits lims) {
    // Our grid element intersects the region limits
    if(lims.x0 < xc && lims.x1 > xc || lims.y0 < yc && lims.y1 > yc) {

      // It intersects all limits (the 4-intersection exists in the element)
      if(lims.x0 < xc && lims.x1 > xc && lims.y0 < yc && lims.y1 > yc) {
        return (vals.bot_left * (xc - lims.x0) * (yc - lims.y0) +
                 vals.top_left * (xc - lims.x0) * (lims.y1 - yc) +
                 vals.top_right * (lims.x1 - xc) * (lims.y1 - yc) +
                 vals.bot_right * (lims.x1 - xc) * (yc - lims.y0)) /
               ((lims.x1 - lims.x0) * (lims.y1 - lims.y0));
      }

      // Intersects the x-limit
      if(lims.x0 < xc && lims.x1 > xc) {

        // It's on the bottom
        if(lims.y0 < xc)
          return (vals.bot_left * (xc - lims.x0) +
                   vals.bot_right * (lims.x1 - xc)) /
                 (lims.x1 - lims.x0);
        // It's on the top
        else
          return (vals.top_left * (xc - lims.x0) +
                   vals.top_right * (lims.x1 - xc)) /
                 (lims.x1 - lims.x0);
      }
      // Intersects the y-limit
      else if(lims.y0 < yc && lims.y1 > yc) {

        // It's on the left
        if(lims.x0 < xc)
          return (vals.bot_left * (yc - lims.y0) +
                   vals.top_left * (lims.y1 - yc)) /
                 (lims.y1 - lims.y0);

        // It's on the right
        else
          return (vals.bot_right * (yc - lims.y0) +
                   vals.top_right * (lims.y1 - yc)) /
                 (lims.y1 - lims.y0);
      }

      // Impossible
      else {
        flog_fatal("Impossible region");
      }
    }

    // Our grid element exists purely in one region
    else {
      if((lims.x0 <= xc) && (lims.y0 <= yc))
        return vals.bot_left;
      if((lims.x0 <= xc) && (lims.y0 >= yc))
        return vals.top_left;
      if((lims.x0 >= xc) && (lims.y0 >= yc))
        return vals.top_right;
      if((lims.x0 >= xc) && (lims.y0 <= yc))
        return vals.bot_right;

      flog_fatal("Region is out of limits");
    }
  };

  // Initial conditions for the problem
  if constexpr(Dim == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {

        const limits lims{
          .x0 = m.template right<ax::x>(i),
          .x1 = m.template left<ax::x>(i),
          .y0 = m.template right<ax::y>(j),
          .y1 = m.template left<ax::y>(j),
        };

        const double ux{find_vals(ux_vals, lims)};
        const double uy{find_vals(uy_vals, lims)};

        mass_density(i, j) = find_vals(rho_vals, lims);
        momentum_density(i, j).x() = mass_density(i, j) * ux;
        momentum_density(i, j).y() = mass_density(i, j) * uy;

        total_energy_density(i, j) =
          mass_density(i, j) *
          (util::find_sie(eos, mass_density(i, j), find_vals(p_vals, lims)) +
            0.5 * (utils::sqr(ux) + utils::sqr(uy)));

      } // for
    }; // forall
  }
  else {
    flog_fatal("Four shocks jet is only defined for D == 2")
  };
} // four_shock_jet

} // namespace hard::tasks::initial_data

#endif // HARD_MODULES_HYDRO_FOUR_SHOCK_JET_HH
