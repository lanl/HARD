#ifndef HARD_APPS_ACTIONS_INITIALIZE_MESH_HH
#define HARD_APPS_ACTIONS_INITIALIZE_MESH_HH

#include "spec/tasks/boundaries/boundary.hh"
#include "spec/types.hh"

namespace hard::actions {

struct init_timestep_region : annotation_action_region {
  inline static const std::string name{"init_timestep"};
};

struct init_boundaries_region : annotation_action_region {
  inline static const std::string name{"init_boundaries"};
};

struct init_topologies_region : annotation_action_region {
  inline static const std::string name{"init_topologies"};
};

struct init_mesh_region : annotation_action_region {
  inline static const std::string name{"init_mesh"};
};

struct init_reconstructor_region : annotation_action_region {
  inline static const std::string name{"init_reconstructor"};
};

struct init_solver_region : annotation_action_region {
  inline static const std::string name{"init_solver"};
};

template<std::size_t D>
void
init_reconstructor(state<D> & s, spec::config_py & config) {
  /*--------------------------------------------------------------------------*
    Reconstruction scheme configuration.
    *--------------------------------------------------------------------------*/
  annotation_rguard<init_reconstructor_region> guard;

  // Default reconstructor is weno5z (set in state initialization)
  if(config.contains("reconstructor")) {
    auto scheme = config["reconstructor"].cast<std::string>();
    if(scheme == "weno5z") {
      s.reconstructor = spec::rc::weno5z;
    }
    else if(scheme == "ppm4") {
      s.reconstructor = spec::rc::ppm4;
    }
    else if(scheme == "minmod") {
      s.reconstructor = spec::rc::minmod;
    }
    else if(scheme == "genminmod") {
      s.reconstructor = spec::rc::genminmod;
    }
    else if(scheme == "noop") {
      s.reconstructor = spec::rc::noop;
    }
    else {
      flog_fatal("unsupported reconstructor(" << scheme << ")");
    }
  }
}

template<std::size_t D>
void
init_timestep(state<D> & s, flecsi::scheduler & sc, spec::config_py & config) {

  annotation_rguard<init_timestep_region> guard;

  /*--------------------------------------------------------------------------*
     Initialize time advance.
    *--------------------------------------------------------------------------*/

  annotation_rguard<tasks::hydro::conserved_to_primitive_task>(),
    sc.execute<tasks::hydro::conserved_to_primitive<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.prim.velocity(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.prim.sound_speed(*s.m),
      s.eos);

  auto lmax_f =
    (annotation_rguard<tasks::hydro::update_max_characteristic_speed_task>(),
      sc.execute<tasks::hydro::update_max_characteristic_speed<D>>(
        flecsi::exec::on,
        *s.m,
        s.cons.hydro.mass_density(*s.m),
        s.prim.velocity(*s.m),
        s.prim.sound_speed(*s.m)));

  s.dtmin_ = (annotation_rguard<tasks::hydro::update_dtmin_task>(),
    sc.reduce<tasks::hydro::update_dtmin<D>, flecsi::exec::fold::min>(
      flecsi::exec::on, *s.m, lmax_f));

  apply_tagged_boundaries<D>(sc, s, tag::boundary);

  /*--------------------------------------------------------------------------*
    Initialize time to 0
   *--------------------------------------------------------------------------*/
  annotation_rguard<tasks::init::init_time_task>(),
    sc.execute<tasks::init::init_time>(
      flecsi::exec::on, s.t(*s.gt), config["t0"].cast<double>());
}

template<std::size_t D>
void
init_solver(state<D> & s, spec::config_py & config) {
  /*--------------------------------------------------------------------------*
    ODE solver scheme configuration.
    *--------------------------------------------------------------------------*/
  annotation_rguard<init_solver_region> guard;

  // Default reconstructor is heun (set in state initialization)
  if(config.contains("ode_solver")) {
    auto scheme = config["ode_solver"].cast<std::string>();
    if(scheme == "rk4") {
      s.use_solver = spec::slv::rk4;
    }
    else if(scheme == "heun") {
      s.use_solver = spec::slv::heun;
    }
    else if(scheme == "euler") {
      s.use_solver = spec::slv::euler;
    }
    else if(scheme == "rk23") {
      s.use_solver = spec::slv::rk23;
    }
    else if(scheme == "rk45") {
      s.use_solver = spec::slv::rk45;
    }
    else {
      flog_fatal("unsupported solver(" << scheme << ")");
    }
  }
}

template<std::size_t D>
typename mesh<D>::periodic_axes
init_boundaries(state<D> & s, spec::config_py & config) {
  annotation_rguard<init_boundaries_region> guard;
  s.bcs = spec::parse_boundaries<D>(config, s.fields());
  return s.bcs.periodic;
}

template<std::size_t D>
void
init_topologies(state<D> & s, flecsi::scheduler & sc) {
  annotation_rguard<init_topologies_region> guard;

  /*--------------------------------------------------------------------------*
    Global and color topology allocations.
   *--------------------------------------------------------------------------*/
  const auto num_colors =
    opt::colors.value() == 0 ? sc.runtime().processes() : opt::colors.value();
  sc.allocate(s.gt, num_colors);
  sc.allocate(s.ct, {num_colors});
}

template<std::size_t D>
void
init_mesh(state<D> & s,
  flecsi::scheduler & sc,
  typename mesh<D>::periodic_axes const & bf,
  spec::config_py & config) {
  using namespace common::utils;

  annotation_rguard<init_mesh_region> guard;

  /*--------------------------------------------------------------------------*
    Mesh topology allocation.
   *--------------------------------------------------------------------------*/

  // Find out how many levels we can have.
  auto get_resolution = [&config](const int dim) {
    return opt::resolution.value() == 0
             ? config["levels"][dim].cast<std::size_t>()
             : opt::resolution.value();
  };

  // Find highest level
  s.min_highest_level = get_resolution(0);
  if(D == 2 || D == 3) {
    s.min_highest_level = std::min(get_resolution(1), s.min_highest_level);
  } // if
  if(D == 3) {
    s.min_highest_level = std::min(get_resolution(2), s.min_highest_level);
  } // if

  // Record lowest level
  if(opt::resolution.value() != 0)
    s.lowest_level = opt::resolution.value();
  else if(config.contains("lowest_level"))
    s.lowest_level = config["lowest_level"].cast<std::size_t>();
  else
    s.lowest_level = s.min_highest_level;

  s.max_num_levels = s.min_highest_level - s.lowest_level + 1;

  if(s.lowest_level > s.min_highest_level)
    flog_fatal("Error in levels setup: lowest_level("
               << s.lowest_level << ") > (" << s.min_highest_level << ")"
               << '\n');

  std::optional<color_distribution> cd;
  if(config.contains("color_distribution")) {
    cd = [&sc,
           cdcfg = config["color_distribution"].cast<color_distribution>()]() {
      return (FLECSI_BACKEND == FLECSI_BACKEND_legion) ||

                 axes_colors<D>(cdcfg) == sc.runtime().processes()
               ? std::optional<color_distribution>(cdcfg)
               : std::nullopt;
    }();
  } // if

  {
    typename mesh<D>::grect geom;
    geom[0][0] = config["coords"][0][0].cast<double>();
    geom[0][1] = config["coords"][1][0].cast<double>();
    if(D == 2 || D == 3) {
      geom[1][0] = config["coords"][0][1].cast<double>();
      geom[1][1] = config["coords"][1][1].cast<double>();
    } // if
    if(D == 3) {
      geom[2][0] = config["coords"][0][2].cast<double>();
      geom[2][1] = config["coords"][1][2].cast<double>();
    } // if

    /*-------------------------------------------------------------------------*
      Set mesh resolution.
     *------------------------------------------------------------------------*/

    for(std::size_t i{0}; i < s.max_num_levels; i++) {
      typename mesh<D>::gcoord axis_extents(D);
      axis_extents[ax::x] = 1 << (get_resolution(0) - i);
      if(D == 2 || D == 3) {
        axis_extents[ax::y] = 1 << (get_resolution(1) - i);
      } // if
      if(D == 3) {
        axis_extents[ax::z] = 1 << (get_resolution(2) - i);
      } // if

      // Add a new grid - the finest grid is already there
      if(i > 0) {
        s.mh.emplace_back(typename mesh<D>::ptr());
      } // if

      if(cd.has_value()) {
        sc.allocate(
          s.mh[i], mesh<D>::color(cd.value(), axis_extents, bf), geom);
      }
      else {
        sc.allocate(s.mh[i],
          mesh<D>::color(sc.runtime().processes(), axis_extents, bf),
          geom);
      }
    } // for
  } // scope
}

} // namespace hard::actions

#endif // HARD_APPS_ACTIONS_INITIALIZE_MESH_HH
