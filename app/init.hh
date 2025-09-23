#ifndef HARD_INITIALIZE_HH
#define HARD_INITIALIZE_HH

#include "options.hh"
#include "spec/eos.hh"
#include "state.hh"
#include "tasks/boundaries/boundary.hh"
#include "tasks/hydro/cons2prim.hh"
#include "tasks/hydro/maxcharspeed.hh"
#include "tasks/init.hh"
#include "tasks/initial_data/all_initial_data.hh"
#include "tasks/io.hh"
#include "tasks/rad.hh"
#include "types.hh"
#include "utils.hh"

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
  flecsi::scheduler & sc = cp.scheduler();

  YAML::Node config = YAML::LoadFile(opt::config.value());

  /*--------------------------------------------------------------------------*
    Solver.
   *--------------------------------------------------------------------------*/

#if defined(USE_FLECSOLVE) && USE_FLECSOLVE
  s.solver_settings.maxiter = config["linear_solver"]["maxiter"].IsDefined()
                                ? config["linear_solver"]["maxiter"].as<int>()
                                : 50;
  s.solver_settings.rtol = config["linear_solver"]["rtol"].IsDefined()
                             ? config["linear_solver"]["rtol"].as<double>()
                             : 1e-12;
  s.solver_settings.use_zero_guess =
    config["linear_solver"]["use_zero_guess"].IsDefined()
      ? config["linear_solver"]["use_zero_guess"].as<bool>()
      : true;
  s.flecsolve_coarse_grid =
    config["linear_solver"]["flecsolve_coarse_grid"].IsDefined()
      ? config["linear_solver"]["flecsolve_coarse_grid"].as<bool>()
      : true;
#endif
  s.jacobi_iterations =
    config["linear_solver"]["jacobi_iterations"].IsDefined()
      ? config["linear_solver"]["jacobi_iterations"].as<double>()
      : 100;

  s.full_multigrid = config["linear_solver"]["full_multigrid"].IsDefined()
                       ? config["linear_solver"]["full_multigrid"].as<bool>()
                       : true;

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

  sc.allocate(s.dense_topology, n_lines);

  const auto num_colors =
    opt::colors.value() == 0 ? sc.runtime().processes() : opt::colors.value();
  sc.allocate(s.gt, num_colors);
  sc.allocate(s.ct, {num_colors});

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

  auto bf =
    execute<tasks::init_boundaries<D>>(flecsi::exec::on, s.bmap(*s.gt), bnds);

  /*--------------------------------------------------------------------------*
    T boundary.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::set_t_boundary>(
    flecsi::exec::on, s.time_boundary(*s.dense_topology), time);
  execute<tasks::init::set_t_boundary>(
    flecsi::exec::on, s.temperature_boundary(*s.dense_topology), temperature);
  if(config["problem"].as<std::string>() == "implosion")
    execute<tasks::init::convert_temperature>(flecsi::exec::on,
      s.temperature_boundary(*s.dense_topology),
      config["temperature_units"].as<std::string>());

    /*--------------------------------------------------------------------------*
      Kappa.
     *--------------------------------------------------------------------------*/
#ifdef ENABLE_RADIATION
  execute<tasks::init::kappa>(s.kappa(*s.gt), config["kappa"].as<double>());
#endif

  /*--------------------------------------------------------------------------*
    Adaptive FLD Check, Closure ID and Limiter ID
   *--------------------------------------------------------------------------*/

#ifdef ENABLE_RADIATION
  // Default is limiter = 1 and closure = 3
  std::size_t ci = config["closure_id"].IsDefined()
                     ? config["closure_id"].as<std::size_t>()
                     : 3;
  std::size_t li = config["limiter_id"].IsDefined()
                     ? config["limiter_id"].as<std::size_t>()
                     : 1;
  sc.execute<tasks::init::closure_id>(s.closure_id(*s.gt), ci);
  sc.execute<tasks::init::limiter_id>(s.limiter_id(*s.gt), li);
#endif

  /*--------------------------------------------------------------------------*
    Gravity Acceleration
   *--------------------------------------------------------------------------*/

  vec<D> g(0.0);
  if(config["gravity_acc"].IsDefined()) {
    g[0] = config["gravity_acc"][0].as<double>();
    if constexpr(D > 1)
      g[1] = config["gravity_acc"][1].as<double>();
    if constexpr(D > 2)
      g[2] = config["gravity_acc"][2].as<double>();
  }
  execute<tasks::init::intialize_gravity_acc<D>>(s.gravity_acc(*s.gt), g);

  /*--------------------------------------------------------------------------*
    Particle mass
   *--------------------------------------------------------------------------*/
  execute<tasks::init::particle_mass>(
    s.particle_mass(*s.gt), config["mean_molecular_weight"].as<double>());

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

  std::optional<color_distribution> cd;
  if(config["color_distribution"]) {
    cd = [&sc,
           cdcfg = config["color_distribution"].as<color_distribution>()]() {
      return (FLECSI_BACKEND == FLECSI_BACKEND_legion) ||

                 util::axes_colors<D>(cdcfg) == sc.runtime().processes()
               ? std::optional<color_distribution>(cdcfg)
               : std::nullopt;
    }();
  } // if

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
        sc.allocate(s.mh[i],
          typename mesh<D>::mpi_coloring(
            sc, cd.value(), axis_extents, bf.get()),
          geom);
      }
      else {
        sc.allocate(s.mh[i],
          typename mesh<D>::mpi_coloring(
            sc, sc.runtime().processes(), axis_extents, bf.get()),
          geom);
      }
    } // for
  } // scope

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

  execute<tasks::init::initialize_gravity_force<D>>(
    flecsi::exec::on, s.gravity_force(*s.m));

  if(config["problem"].as<std::string>() == "sod") {

#ifdef ENABLE_RADIATION
    flog_fatal("Sod must be built with ENABLE_RADIATION=OFF");
#endif

    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::sod, D>>(
      flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "rankine-hugoniot") {

#ifdef ENABLE_RADIATION
    flog_fatal("Rankine-Hugoniot must be built with ENABLE_RADIATION=OFF");
#endif

    execute<tasks::initial_data::
        shock<tasks::initial_data::shock_tubes::rankine_hugoniot, D>>(
      flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "leblanc") {

#ifdef ENABLE_RADIATION
    flog_fatal("Leblanc must be built with ENABLE_RADIATION=OFF");
#endif

    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::leblanc, D>>(
      flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "acoustic-wave") {

#ifdef ENABLE_RADIATION
    flog_fatal("Acoustic wave must be built with ENABLE_RADIATION=OFF");
#endif

    execute<tasks::initial_data::acoustic_wave<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "kh-test") {

#ifdef ENABLE_RADIATION
    flog_fatal(
      "Kelvin-Helmholtz instability must be built with ENABLE_RADIATION=OFF");
#endif

    execute<tasks::initial_data::kh_instability<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  // Rayleigh-Taylor setup
  else if(config["problem"].as<std::string>() == "rt-test") {

#ifdef ENABLE_RADIATION
    flog_fatal(
      "Rayleigh-Taylor instability must be built with ENABLE_RADIATION=OFF");
#endif

    execute<tasks::initial_data::rt_instability<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.gravity_force(*s.m),
      s.gravity_acc(*s.gt),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "heating_and_cooling") {
    if(config["eos"].as<std::string>() != "ideal")
      flog_fatal("Heating and cooling test only supports Ideal Gas eos");
    execute<tasks::initial_data::heating_and_cooling<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.particle_mass(*s.gt),
      config["gamma"].as<double>());
  }
  // Heating and Cooling for AFLD
  else if(config["problem"].as<std::string>() == "heating-cooling-afld") {
    if(config["eos"].as<std::string>() != "ideal")
      flog_fatal("Heating and cooling test only supports Ideal Gas eos");
    sc.execute<tasks::initial_data::heating_and_cooling_afld<D>>(
      flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.particle_mass(*s.gt),
      config["gamma"].as<double>());
  }
  else if(config["problem"].as<std::string>() == "sedov") {

#ifdef ENABLE_RADIATION
    flog_fatal("Sedov blast must be built with ENABLE_RADIATION=OFF");
#endif

    execute<tasks::initial_data::sedov_blast<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m));
  }
  else if(config["problem"].as<std::string>() == "implosion") {
    execute<tasks::initial_data::implosion_forced_T<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.temperature_boundary(*s.dense_topology),
      s.particle_mass(*s.gt),
      config["gamma"].as<double>());
    s.mg = true;
  }
  // FIXME: This problem has not been tested for correctness
  else if(config["problem"].as<std::string>() == "rad-rh") {
    execute<tasks::initial_data::
        rad_RH<tasks::initial_data::rad_shock::rad_rankine_hugoniot, D>>(
      flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      config["gamma"].as<double>(),
      s.particle_mass(*s.gt));
  }
  // FIXME: This problem has not been tested for correctness
  else if(config["problem"].as<std::string>() == "lw-implosion") {

#ifdef ENABLE_RADIATION
    flog_fatal(
      "Liska-Wendroff implosion must be built with ENABLE_RADIATION=OFF");
#endif

    execute<tasks::initial_data::lw_implosion<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      config["gamma"].as<double>());
  }
  // Kelvin Helmholtz with radiation setup
  else if(config["problem"].as<std::string>() == "kh-rad-test") {

    execute<tasks::initial_data::kh_instability_rad<D>>(flecsi::exec::on,
      *s.m,
      s.mass_density(*s.m),
      s.momentum_density(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m),
      s.eos);
  }
  else {
    flog_fatal(
      "unsupported problem(" << config["problem"].as<std::string>() << ")");
  } // if

  /*--------------------------------------------------------------------------*
    Initialize time advance.
   *--------------------------------------------------------------------------*/

  sc.execute<tasks::hydro::conservative_to_primitive<D>>(flecsi::exec::on,
    *s.m,
    s.mass_density(*s.m),
    s.momentum_density(*s.m),
    s.total_energy_density(*s.m),
    s.velocity(*s.m),
    s.pressure(*s.m),
    s.specific_internal_energy(*s.m),
    s.sound_speed(*s.m),
    s.eos);
  auto lmax_f = sc.execute<tasks::hydro::update_max_characteristic_speed<D>>(
    flecsi::exec::on,
    *s.m,
    s.mass_density(*s.m),
    s.velocity(*s.m),
    s.sound_speed(*s.m));
  s.dtmin_ = reduce<hard::task::rad::update_dtmin<D>, exec::fold::min>(
    flecsi::exec::on, *s.m, lmax_f);

  sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
    *s.m,
    s.bmap(*s.gt),
    std::vector{s.mass_density(*s.m),
      s.pressure(*s.m),
      s.specific_internal_energy(*s.m),
      s.radiation_energy_density(*s.m),
      s.total_energy_density(*s.m)},
    std::vector{s.velocity(*s.m), s.momentum_density(*s.m)});
  if(s.mg) {
    // FIXME: figure out how not to use the hardcoded radiation temperature
    // boundary
    auto radiation_boundary_f =
      sc.execute<task::rad::interp_e_boundary>(flecsi::exec::on,
        s.t(*s.gt),
        s.time_boundary(*s.dense_topology),
        s.temperature_boundary(*s.dense_topology));
    sc.execute<tasks::apply_dirichlet_boundaries<D>>(flecsi::exec::on,
      *s.m,
      s.bmap(*s.gt),
      std::vector{s.radiation_energy_density(*s.m)},
      radiation_boundary_f);
  }

  /*--------------------------------------------------------------------------*
    Initialize time to 0
   *--------------------------------------------------------------------------*/
  sc.execute<tasks::init::init_time>(
    flecsi::exec::on, s.t(*s.gt), config["t0"].as<double>());

  /*--------------------------------------------------------------------------*
    Save raw data after initialization, at t=t_i
  *--------------------------------------------------------------------------*/
#ifndef HARD_BENCHMARK_MODE

  auto lm = data::launch::make(sc, *s.m);
  execute<tasks::io::csv<D>, mpi>(flecsi::exec::on,
    spec::io::name{""} << std::setfill('0') << std::setw(5) << cp.step(),
    s.t(*s.gt),
    lm,
    std::vector{s.mass_density(lm),
      s.pressure(lm),
      s.sound_speed(lm),
      s.specific_internal_energy(lm),
      s.total_energy_density(lm),
      s.radiation_energy_density(lm)},
    std::vector{s.velocity(lm), s.momentum_density(lm)},
    std::vector<std::string>{"density",
      "pressure",
      "sound_speed",
      "specific_internal_energy",
      "total_energy_density",
      "radiation_energy_density"},
    std::vector<std::string>{"velocity", "momentum_density"});
#endif
  /*--------------------------------------------------------------------------*
    Initialize catalyst: Build the lattice, initialize adaptor, parse yaml file
    Then send initial problem state down the conduit
   *--------------------------------------------------------------------------*/

#ifdef USE_CATALYST
  if constexpr(D == 3) {
    // Initialize catalyst
    execute<tasks::external::initialize, mpi>(flecsi::exec::on);
    sc.allocate(s.pt, num_colors);

    // Initialize lattice
    execute<tasks::external::get_lattice_size<D>, mpi>(flecsi::exec::on, *s.m);
    flog(info) << "init action, catalyst: got lattice dimensions from hard: "
               << number_vertices[0] << " x " << number_vertices[1] << " x "
               << number_vertices[2] << " vertices." << std::endl;

    execute<tasks::external::get_color_data<D>, mpi>(flecsi::exec::on, *s.m);
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
    execute<tasks::external::init_attributes, mpi>(flecsi::exec::on,
      s.catalyst_data(*s.pt),
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
    execute<tasks::external::update_attributes<D>, mpi>(flecsi::exec::on,
      s.catalyst_data(*s.pt),
      s.t(*s.gt),
      *s.m,
      s.mass_density(*s.m),
      s.velocity(*s.m),
      s.pressure(*s.m),
      s.total_energy_density(*s.m),
      s.radiation_energy_density(*s.m)); // <<< add variables here for catalyst
    flog(info) << "init action, catalyst: execute catalyst for initial state"
               << std::endl;
    execute<tasks::external::execute_catalyst, mpi>(
      flecsi::exec::on, s.catalyst_data(*s.pt), 0, s.t(*s.gt), lattice);
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
