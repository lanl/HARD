#ifndef FLASTRO_INITIALIZE_HH
#define FLASTRO_INITIALIZE_HH

#include "options.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/hydro.hh"
#include "tasks/hydro/cons2prim.hh"
#include "tasks/init.hh"
#include "tasks/initial_data/all_initial_data.hh"
#include "tasks/utils.hh"
#include "types.hh"

#ifdef USE_CATALYST
#include "tasks/catalyst.hh"
#endif

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

namespace flastro {
namespace action {

template<std::size_t D>
void
initialize(control_policy<state, D> & cp) {
  using namespace flecsi;
  auto & s = cp.state();
  YAML::Node config = YAML::LoadFile(opt::config.value());

  /*--------------------------------------------------------------------------*
    Global and color topology allocations.
   *--------------------------------------------------------------------------*/

  s.gt.allocate({});
  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();
  s.ct.allocate(num_colors);

  /*--------------------------------------------------------------------------*
    Set boundaries.
   *--------------------------------------------------------------------------*/

  std::array<std::array<bd::boundary_type, 2>, D> bnds;
  bnds[ax::x][bd::low] =
    utils::mesh_boundary<D>(config["boundaries"]["xlow"].as<std::string>());
  bnds[ax::x][bd::high] =
    utils::mesh_boundary<D>(config["boundaries"]["xhigh"].as<std::string>());
  if(D == 2 || D == 3) {
    bnds[ax::y][bd::low] =
      utils::mesh_boundary<D>(config["boundaries"]["ylow"].as<std::string>());
    bnds[ax::y][bd::high] =
      utils::mesh_boundary<D>(config["boundaries"]["yhigh"].as<std::string>());
  } // if
  if(D == 3) {
    bnds[ax::z][bd::low] =
      utils::mesh_boundary<D>(config["boundaries"]["zlow"].as<std::string>());
    bnds[ax::z][bd::high] =
      utils::mesh_boundary<D>(config["boundaries"]["zhigh"].as<std::string>());
  } // if

  auto bf = execute<tasks::init::boundaries<D>>(s.bmap(s.gt), bnds);

  /*--------------------------------------------------------------------------*
    Gamma.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::gamma>(gamma(s.gt), config["gamma"].as<double>());

  /*--------------------------------------------------------------------------*
    Kappa.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::kappa>(kappa(s.gt), config["kappa"].as<double>());

  /*--------------------------------------------------------------------------*
    Particle mass
   *--------------------------------------------------------------------------*/
  execute<tasks::init::particle_mass>(
    particle_mass(s.gt), config["mean_molecular_weight"].as<double>());

  /*--------------------------------------------------------------------------*
    Mesh topology allocation.
   *--------------------------------------------------------------------------*/

  // Find out how many levels we can have.

  // Record lowest level
  s.lowest_level = config["lowest_level"].as<std::size_t>();

  // Find highest level
  s.highest_level = config["levels"][0].as<std::size_t>();
  if(D == 2 || D == 3) {
    s.highest_level =
      std::max(s.highest_level, config["levels"][1].as<std::size_t>());
  } // if
  if(D == 3) {
    s.highest_level =
      std::max(s.highest_level, config["levels"][2].as<std::size_t>());
  } // if
  s.max_num_levels = s.highest_level - s.lowest_level + 1;

  {
    typename mesh<D>::grect geom;
    geom[0][0] = config["coords"][0][0].as<double>();
    geom[0][1] = config["coords"][1][0].as<double>();
    if(D == 2 || D == 3) {
      geom[1][0] = config["coords"][0][1].as<double>();
      geom[1][1] = config["coords"][1][1].as<double>();
    } // if
    if(D == 3) {
      geom[2][0] = config["coords"][0][2].as<double>();
      geom[2][1] = config["coords"][1][2].as<double>();
    } // if

    /*-------------------------------------------------------------------------*
      Set mesh resolution.
     *------------------------------------------------------------------------*/

    for(std::size_t i{0}; i < s.max_num_levels; i++) {
      typename mesh<D>::gcoord axis_extents(D);
      axis_extents[ax::x] =
        std::pow(2, config["levels"][0].as<std::size_t>() - i);
      if(D == 2 || D == 3) {
        axis_extents[ax::y] =
          std::pow(2, config["levels"][1].as<std::size_t>() - i);
      } // if
      if(D == 3) {
        axis_extents[ax::z] =
          std::pow(2, config["levels"][2].as<std::size_t>() - i);
      } // if

      // Add a new grid - the finest grid is already there
      if(i > 0) {
        s.mh.emplace_back(std::make_unique<typename mesh<D>::slot>());
      } // if

      // Allocate the grid
      s.mh[i]->allocate(
        typename mesh<D>::mpi_coloring{num_colors, axis_extents, bf.get()},
        geom);
    } // for
  } // scope

  /*--------------------------------------------------------------------------*
    Fake initialization to avoid legion warnings.
   *--------------------------------------------------------------------------*/

  // clang-format off
  execute<tasks::init::touch<D>>(s.m,
    s.mass_density(s.m), s.momentum_density(s.m), s.total_energy_density(s.m),
    s.velocity(s.m), s.pressure(s.m),
    s.rTail(s.m), s.ruTail(s.m), s.rETail(s.m), s.uTail(s.m), s.pTail(s.m), 
    s.rHead(s.m), s.ruHead(s.m), s.rEHead(s.m), s.uHead(s.m), s.pHead(s.m),
    s.rF(s.m), s.ruF(s.m), s.rEF(s.m)
#ifndef DISABLE_RADIATION
    , s.radiation_energy_density(s.m), s.EradTail(s.m), s.EradHead(s.m), s.EradF(s.m),
    s.Esf(s.m), s.Ef(s.m), s.Ew(s.m), s.Diff(s.m), s.Resf(s.m), s.Errf(s.m),
    s.gradient_rad_energy(s.m), s.magnitude_gradient_rad_energy(s.m),
    s.radiation_force(s.m), s.R_value(s.m), s.lambda_bridge(s.m),
    s.velocity_gradient(s.m)
#endif
    );
  execute<tasks::init::touch1<D>>(s.m,
    s.mass_density(s.m), s.momentum_density(s.m), s.total_energy_density(s.m),
    s.velocity(s.m), s.pressure(s.m),
    s.rTail(s.m), s.ruTail(s.m), s.rETail(s.m), s.uTail(s.m), s.pTail(s.m), 
    s.rHead(s.m), s.ruHead(s.m), s.rEHead(s.m), s.uHead(s.m), s.pHead(s.m),
    s.rF(s.m), s.ruF(s.m), s.rEF(s.m)
#ifndef DISABLE_RADIATION
    , s.radiation_energy_density(s.m), s.EradTail(s.m), s.EradHead(s.m), s.EradF(s.m),
    s.Esf(s.m), s.Ef(s.m), s.Ew(s.m), s.Diff(s.m), s.Resf(s.m), s.Errf(s.m),
    s.gradient_rad_energy(s.m), s.magnitude_gradient_rad_energy(s.m),
    s.radiation_force(s.m), s.R_value(s.m), s.lambda_bridge(s.m),
    s.velocity_gradient(s.m)
#endif
    );
  // clang-format on

  /*--------------------------------------------------------------------------*
    Initialize problem state.
   *--------------------------------------------------------------------------*/

  if(config["problem"].as<std::string>() == "sod") {
    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::sod, D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt));
  }
  else if(config["problem"].as<std::string>() == "rankine-hugoniot") {
    execute<tasks::initial_data::
              shock<tasks::initial_data::shock_tubes::rankine_hugoniot, D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt));
  }
  else if(config["problem"].as<std::string>() == "sine-wave") {
    execute<tasks::initial_data::sine_wave<D>, flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt));
  }
  else if(config["problem"].as<std::string>() == "kh-test") {
    execute<tasks::initial_data::kh_instability<D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt));
  }
  else if(config["problem"].as<std::string>() == "radiative-kh-instability") {
    execute<tasks::initial_data::radiative_kh_instability<D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      particle_mass(s.gt));
  }
  else if(config["problem"].as<std::string>() == "heating_and_cooling") {
    execute<tasks::initial_data::heating_and_cooling<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      particle_mass(s.gt));
  }
  else if(config["problem"].as<std::string>() == "rad-rh") {
    execute<tasks::initial_data::rad_RH<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      particle_mass(s.gt));
  }
  else if(config["problem"].as<std::string>() == "sedov") {
    execute<tasks::initial_data::sedov_blast<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      particle_mass(s.gt));
  }
  else if(config["problem"].as<std::string>() == "shadow_2d") {
    execute<tasks::initial_data::shadow_2d<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt));
  }
  else if(config["problem"].as<std::string>() == "radiation_energy_diffusion") {
    execute<tasks::initial_data::radiation_energy_diffusion<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m));
  }
  else if(config["problem"].as<std::string>() == "rayleigh_taylor") {
    execute<tasks::initial_data::rayleigh_taylor<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      particle_mass(s.gt));
  }
  else if(config["problem"].as<std::string>() == "su-olson") {
    execute<tasks::initial_data::su_olson<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m));
  }
  else {
    flog_fatal(
      "unsupported problem(" << config["problem"].as<std::string>() << ")");
  } // if

  /*--------------------------------------------------------------------------*
    Initialize time advance.
   *--------------------------------------------------------------------------*/

  execute<tasks::apply_boundaries<D>, flecsi::default_accelerator>(s.m,
    s.bmap(s.gt),
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.radiation_energy_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m));
  execute<tasks::hydro::conservative_to_primitive<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    gamma(s.gt));
  auto lmax_f = execute<tasks::hydro::update_max_characteristic_speed<D>,
    flecsi::default_accelerator>(
    s.m, s.mass_density(s.m), s.velocity(s.m), s.pressure(s.m), gamma(s.gt));
  s.dtmin_ =
    reduce<tasks::hydro::update_dtmin<D>, exec::fold::min>(s.m, lmax_f);

  /*--------------------------------------------------------------------------*
    Initialize time to 0
   *--------------------------------------------------------------------------*/
  flecsi::execute<tasks::init::init_time>(s.t(s.gt), config["t0"].as<double>());

  /*--------------------------------------------------------------------------*
    Save raw data after initialization, at t=t_i
  *--------------------------------------------------------------------------*/
#ifndef FLASTRO_BENCHMARK_MODE

  auto lm = data::launch::make(s.m);
  execute<tasks::io::raw<D>, mpi>(
    spec::io::name{"output-"} << std::setfill('0') << std::setw(5) << cp.step(),
    s.t(s.gt),
    lm,
    s.mass_density(lm),
    s.pressure(lm),
    s.velocity(lm),
    s.momentum_density(lm),
    s.total_energy_density(lm),
    s.radiation_energy_density(lm));

#endif
  /*--------------------------------------------------------------------------*
    Initialize catalyst: Build the lattice, initialize adaptor, parse yaml file
    Then send initial problem state down the conduit
   *--------------------------------------------------------------------------*/

#ifdef USE_CATALYST
  if constexpr(D == 3) {
    // Initialize catalyst
    execute<tasks::external::initialize, mpi>();

    // Initialize lattice
    execute<tasks::external::get_lattice_size<D>, mpi>(s.m);
    flog(info) << "init action, catalyst: got lattice dimensions from flastro: "
               << number_vertices[0] << " x " << number_vertices[1] << " x "
               << number_vertices[2] << " vertices." << std::endl;

    execute<tasks::external::get_color_data<D>, mpi>(s.m);
    flog(info) << "init action, catalyst: got color dimensions from flastro: "
               << number_colors[0] << " x " << number_colors[1] << " x "
               << number_colors[2] << " color blocks." << std::endl;

    flog(info)
      << "init action, catalyst, initialize lattice: number of cells = "
      << number_vertices[0] - 1 << " x " << number_vertices[1] - 1 << " x "
      << number_vertices[2] - 1 << std::endl;

    std::cout << "color: " << my_rank << " x: " << std::setw(5)
              << color_cell_id_begin[0] << " -- " << std::setw(5)
              << color_cell_id_end[0] << " y: " << std::setw(5)
              << color_cell_id_begin[1] << " -- " << std::setw(5)
              << color_cell_id_end[1] << " z: " << std::setw(5)
              << color_cell_id_begin[2] << " -- " << std::setw(5)
              << color_cell_id_end[2] << std::endl;

    flog(info)
      << "init action, catalyst, initialize lattice: lattice spacing = "
      << lattice_spacing[0] << ", " << lattice_spacing[1] << ", "
      << lattice_spacing[2] << std::endl;

    lattice.initialize(number_vertices, lattice_spacing);

    // Check "Cells", just to confirm it was initialized correctly and
    // understand what it is
    flog(info)
      << "init action, catalyst, initialize lattice: this process owns "
      << lattice.get_number_of_cells() << " cells and "
      << lattice.get_number_of_points() << " points (vertices)" << std::endl;
    /*
    std::stringstream info_s;
    info_s << "Points/Vertex coordinates (x,y,z)_1, ..., (x,y,z)_n: ";
    copy(lattice.get_points_array(),lattice.get_points_array() +
    (lattice.get_number_of_points() * 3),std::ostream_iterator<double>(info_s, "
    ")); flog(info) << info_s.str() << std::endl;
    */

    // Initialize catalyst data structure
    execute<tasks::external::init_attributes, mpi>(catalyst_data(pt),
      lattice.get_number_of_points(),
      lattice.get_number_of_cells());

    // Scan yaml file for rendering info
    if(config["catalyst"]["render"]) {
      flog(info) << "init action: catalyst will render "
                 << config["catalyst"]["render"].size()
                 << " fields:" << std::endl;
      for(YAML::const_iterator it = config["catalyst"]["render"].begin();
          it != config["catalyst"]["render"].end();
          ++it) {
        flog(info) << "\t -" << it->as<std::string>() << std::endl;
      }
    }
    else {
      flog(info) << "init action, catalyst: no catalyst rendering information "
                    "provided in yaml file!"
                 << std::endl;
    }

    // Send initial problem state to catalyst
    execute<tasks::external::update_attributes<D>, mpi>(catalyst_data(pt),
      s.t(s.gt),
      s.m,
      s.mass_density(s.m),
      s.velocity(s.m),
      s.pressure(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m)); // <<< add variables here for catalyst
    flog(info) << "init action, catalyst: execute catalyst for initial state"
               << std::endl;
    execute<tasks::external::execute_catalyst, mpi>(
      catalyst_data(pt), 0, s.t(s.gt), lattice);
  }
  else {
    /* Do nothing, catalyst/paraview visualization only for 3D */
  }
#endif

} // initialize

inline control<state, 1>::action<initialize<1>, cp::initialize> init1_action;
inline control<state, 2>::action<initialize<2>, cp::initialize> init2_action;
inline control<state, 3>::action<initialize<3>, cp::initialize> init3_action;

} // namespace action
} // namespace flastro

#endif // FLASTRO_INITIALIZE_HH
