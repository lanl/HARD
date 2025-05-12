#ifndef HARD_INITIALIZE_HH
#define HARD_INITIALIZE_HH

#include "options.hh"
#include "spec/eos.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/hydro/cons2prim.hh"
#include "tasks/hydro/maxcharspeed.hh"
#include "tasks/init.hh"
#include "tasks/initial_data/all_initial_data.hh"
#include "tasks/io.hh"
#include "tasks/rad.hh"
#include "types.hh"

#ifdef USE_CATALYST
#include "tasks/catalyst.hh"
#endif

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

namespace hard {
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

  // How many boundary points to allocate
  std::string filename{opt::source_fds};

  int n_lines{0};
  std::ifstream filein(filename);
  std::vector<double> time;
  std::vector<double> temperature;
  for(std::string line; std::getline(filein, line);) {

    int i{0};
    std::istringstream input;
    input.str(line);

    for(std::string element; std::getline(input, element, ' ');) {
      if(i == 0) {
        time.emplace_back(std::stod(element));
      }
      else {
        temperature.emplace_back(std::stod(element));
      }
      i++;
    }

    n_lines++;
  }

  s.dense_topology.allocate(n_lines);

  s.gt.allocate({});
  const auto num_colors =
    opt::colors.value() == 0 ? flecsi::processes() : opt::colors.value();
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
    T boundary.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::set_t_boundary>(time_boundary(s.dense_topology), time);
  execute<tasks::init::set_t_boundary>(
    temperature_boundary(s.dense_topology), temperature);
  if(config["problem"].as<std::string>() == "implosion")
    execute<tasks::init::convert_temperature>(
      temperature_boundary(s.dense_topology),
      config["temperature_units"].as<std::string>());

    /*--------------------------------------------------------------------------*
      Kappa.
     *--------------------------------------------------------------------------*/
#ifdef ENABLE_RADIATION
  execute<tasks::init::kappa>(kappa(s.gt), config["kappa"].as<double>());
#endif
  /*--------------------------------------------------------------------------*
    Particle mass
   *--------------------------------------------------------------------------*/
  execute<tasks::init::particle_mass>(
    particle_mass(s.gt), config["mean_molecular_weight"].as<double>());

  /*--------------------------------------------------------------------------*
    Mesh topology allocation.
   *--------------------------------------------------------------------------*/

  // Find out how many levels we can have.
  auto get_resolution = [&config](const int dim) {
    return opt::resolution.value() == 0
             ? config["levels"][dim].as<std::size_t>()
             : opt::resolution.value();
  };

  // Record lowest level
  s.lowest_level = opt::resolution.value() == 0
                     ? config["lowest_level"].as<std::size_t>()
                     : opt::resolution.value();

  // Find highest level
  s.highest_level = get_resolution(0);
  if(D == 2 || D == 3) {
    s.highest_level = get_resolution(1);
  } // if
  if(D == 3) {
    s.highest_level = get_resolution(2);
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
      auto opt_dx = opt::resolution.value();
      typename mesh<D>::gcoord axis_extents(D);
      axis_extents[ax::x] = 1 << get_resolution(0) - i;
      if(D == 2 || D == 3) {
        axis_extents[ax::y] = 1 << get_resolution(1) - i;
      } // if
      if(D == 3) {
        axis_extents[ax::z] = 1 << get_resolution(2) - i;
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
#ifdef ENABLE_RADIATION
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
#ifdef ENABLE_RADIATION
    , s.radiation_energy_density(s.m), s.EradTail(s.m), s.EradHead(s.m), s.EradF(s.m),
    s.Esf(s.m), s.Ef(s.m), s.Ew(s.m), s.Diff(s.m), s.Resf(s.m), s.Errf(s.m),
    s.gradient_rad_energy(s.m), s.magnitude_gradient_rad_energy(s.m),
    s.radiation_force(s.m), s.R_value(s.m), s.lambda_bridge(s.m),
    s.velocity_gradient(s.m)
#endif
    );
  // clang-format on

  /*--------------------------------------------------------------------------*
    Equation of State
   *--------------------------------------------------------------------------*/

  if(config["eos"].as<std::string>() == "ideal") {
    s.eos =
      singularity::IdealGas(config["gamma"].as<double>() - 1, 2.0 /* FIXME */);
  }
  else if(config["eos"].as<std::string>() == "spiner") {
    s.eos = singularity::SpinerEOSDependsRhoSie(
      config["spiner_file"].as<std::string>(),
      config["spiner_matid"].as<std::string>());
  }
  else if(config["eos"].as<std::string>() == "gruneisen") {
    s.eos = singularity::Gruneisen(config["gruneisen_c0"].as<double>(),
      config["gruneisen_s1"].as<double>(),
      0., // s2
      0., // s3
      config["gruneisen_G0"].as<double>(),
      0., // b
      0., // rho0
      config["gruneisen_T0"].as<double>(),
      0., // P0
      config["gruneisen_Cv"].as<double>(),
      0. // rho_max

    );
  }
  else {
    flog_fatal("unsupported EOS");
  }

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
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "rankine-hugoniot") {
    execute<tasks::initial_data::
              shock<tasks::initial_data::shock_tubes::rankine_hugoniot, D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "leblanc") {
    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::leblanc, D>,
      flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.eos);
  }

  else if(config["problem"].as<std::string>() == "sine-wave") {
    execute<tasks::initial_data::sine_wave<D>, flecsi::default_accelerator>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "acoustic-wave") {
    execute<tasks::initial_data::acoustic_wave<D>, flecsi::default_accelerator>(
      s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "kh-test") {
    execute<tasks::initial_data::kh_instability<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "heating_and_cooling") {
    if(config["eos"].as<std::string>() != "ideal")
      flog_fatal("Heating and cooling test only supports Ideal Gas eos");
    execute<tasks::initial_data::heating_and_cooling<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      particle_mass(s.gt),
      config["gamma"].as<double>());
  }
  else if(config["problem"].as<std::string>() == "sedov") {
    execute<tasks::initial_data::sedov_blast<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m));
  }
  else if(config["problem"].as<std::string>() == "implosion") {
    execute<tasks::initial_data::implosion_forced_T<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      temperature_boundary(s.dense_topology),
      particle_mass(s.gt),
      config["gamma"].as<double>());
    s.mg = true;
  }
#if 0
  else if(config["problem"].as<std::string>() == "rad-rh") {
    execute<tasks::initial_data::
        rad_RH<tasks::initial_data::rad_shock::rad_rankine_hugoniot, D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt),
      particle_mass(s.gt));
  }
  else if(config["problem"].as<std::string>() == "lw-implosion") {
    execute<tasks::initial_data::lw_implosion<D>>(s.m,
      s.mass_density(s.m),
      s.momentum_density(s.m),
      s.total_energy_density(s.m),
      s.radiation_energy_density(s.m),
      gamma(s.gt));
  }
#endif
  else {
    flog_fatal(
      "unsupported problem(" << config["problem"].as<std::string>() << ")");
  } // if

  /*--------------------------------------------------------------------------*
    Initialize time advance.
   *--------------------------------------------------------------------------*/

  if(s.mg) {
    // FIXME: figure out how not to use the hardcoded radiation temperature
    // boundary
    auto radiation_boundary_f =
      flecsi::execute<task::rad::interp_e_boundary>(s.t(s.gt),
        time_boundary(s.dense_topology),
        temperature_boundary(s.dense_topology));
    flecsi::execute<tasks::apply_radiation_boundary<D>,
      flecsi::default_accelerator>(
      s.m, s.radiation_energy_density(s.m), radiation_boundary_f);
  }
  execute<tasks::hydro::conservative_to_primitive<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.specific_internal_energy(s.m),
    s.sound_speed(s.m),
    s.eos);
  auto lmax_f = execute<tasks::hydro::update_max_characteristic_speed<D>,
    flecsi::default_accelerator>(s.m,
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.sound_speed(s.m));
  s.dtmin_ =
    reduce<hard::task::rad::update_dtmin<D>, exec::fold::min>(s.m, lmax_f);

  execute<tasks::apply_boundaries<D>, flecsi::default_accelerator>(s.m,
    s.bmap(s.gt),
    s.mass_density(s.m),
    s.velocity(s.m),
    s.pressure(s.m),
    s.specific_internal_energy(s.m),
    s.radiation_energy_density(s.m),
    s.momentum_density(s.m),
    s.total_energy_density(s.m));

  /*--------------------------------------------------------------------------*
    Initialize time to 0
   *--------------------------------------------------------------------------*/
  flecsi::execute<tasks::init::init_time>(s.t(s.gt), config["t0"].as<double>());

  /*--------------------------------------------------------------------------*
    Save raw data after initialization, at t=t_i
  *--------------------------------------------------------------------------*/
#ifndef HARD_BENCHMARK_MODE

  auto lm = data::launch::make(s.m);
  execute<tasks::io::raw<D>, mpi>(
    spec::io::name{"output-"} << std::setfill('0') << std::setw(5) << cp.step(),
    s.t(s.gt),
    lm,
    s.mass_density(lm),
    s.pressure(lm),
    s.sound_speed(lm),
    s.specific_internal_energy(lm),
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
    flog(info) << "init action, catalyst: got lattice dimensions from hard: "
               << number_vertices[0] << " x " << number_vertices[1] << " x "
               << number_vertices[2] << " vertices." << std::endl;

    execute<tasks::external::get_color_data<D>, mpi>(s.m);
    flog(info) << "init action, catalyst: got color dimensions from hard: "
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
} // namespace hard

#endif // HARD_INITIALIZE_HH
