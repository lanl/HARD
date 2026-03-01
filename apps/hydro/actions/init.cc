#include "options.hh"
#include "state.hh"

#include "common/utils.hh"
#include "hydro/tasks/cons2prim.hh"
#include "hydro/tasks/init.hh"
#include "hydro/tasks/initial_data/all_initial_data.hh"
#include "hydro/tasks/maxcharspeed.hh"
#include "hydro/tasks/rhs.hh"
#include "spec/eos.hh"
#include "spec/tasks/boundaries/boundary.hh"
#include "spec/tasks/io.hh"

#include <flecsi/flog.hh>
#include <yaml-cpp/yaml.h>

namespace hard {

template<std::size_t D>
void
initialize(control_policy<state, D> & cp) {
  using namespace flecsi;
  using namespace common::utils;
  auto & s = cp.state();
  flecsi::scheduler & sc = cp.scheduler();

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

  auto bf = execute<tasks::init_boundaries<D>>(
    flecsi::exec::on, s.icst.bmap(*s.gt), bnds);

  /*--------------------------------------------------------------------------*
    T boundary.
   *--------------------------------------------------------------------------*/

  execute<tasks::init::set_t_boundary>(
    flecsi::exec::on, s.icst.time_boundary(*s.dense_topology), time);
  execute<tasks::init::set_t_boundary>(flecsi::exec::on,
    s.icst.temperature_boundary(*s.dense_topology),
    temperature);
  if(config["problem"].as<std::string>() == "implosion")
    execute<tasks::init::convert_temperature>(flecsi::exec::on,
      s.icst.temperature_boundary(*s.dense_topology),
      config["temperature_units"].as<std::string>());

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
  execute<tasks::init::initialize_gravity_acc<D>>(s.icst.gravity_acc(*s.gt), g);

  /*--------------------------------------------------------------------------*
    Particle mass
   *--------------------------------------------------------------------------*/
  execute<tasks::init::particle_mass>(
    s.icst.particle_mass(*s.gt), config["mean_molecular_weight"].as<double>());

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
    s.highest_level = std::min(get_resolution(1), s.highest_level);
  } // if
  if(D == 3) {
    s.highest_level = std::min(get_resolution(2), s.highest_level);
  } // if
  s.max_num_levels = s.highest_level - s.lowest_level + 1;

  if(s.lowest_level > s.highest_level)
    flog_fatal("Error in levels setup: lowest_level("
               << s.lowest_level << ") > (" << s.highest_level << ")" << '\n');

  std::optional<color_distribution> cd;
  if(config["color_distribution"]) {
    cd = [&sc,
           cdcfg = config["color_distribution"].as<color_distribution>()]() {
      return (FLECSI_BACKEND == FLECSI_BACKEND_legion) ||

                 axes_colors<D>(cdcfg) == sc.runtime().processes()
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
    flecsi::exec::on, s.src_t.hydro.gravity_force(*s.m));

  if(config["problem"].as<std::string>() == "sod") {
    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::sod, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "rankine-hugoniot") {
    execute<tasks::initial_data::
        shock<tasks::initial_data::shock_tubes::rankine_hugoniot, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "leblanc") {
    execute<
      tasks::initial_data::shock<tasks::initial_data::shock_tubes::leblanc, D>>(
      flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "acoustic-wave") {
    execute<tasks::initial_data::acoustic_wave<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "kh-test") {
    execute<tasks::initial_data::kh_instability<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  // Rayleigh-Taylor setup
  else if(config["problem"].as<std::string>() == "rt-test") {
    execute<tasks::initial_data::rt_instability<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.src_t.hydro.gravity_force(*s.m),
      s.icst.gravity_acc(*s.gt),
      s.cons.hydro.total_energy_density(*s.m),
      s.eos);
  }
  else if(config["problem"].as<std::string>() == "sedov") {
    execute<tasks::initial_data::sedov_blast<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m));
  }
  else if(config["problem"].as<std::string>() == "lw-implosion") {
    execute<tasks::initial_data::lw_implosion<D>>(flecsi::exec::on,
      *s.m,
      s.cons.hydro.mass_density(*s.m),
      s.cons.hydro.momentum_density(*s.m),
      s.cons.hydro.total_energy_density(*s.m),
      config["gamma"].as<double>());
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
    s.cons.hydro.mass_density(*s.m),
    s.cons.hydro.momentum_density(*s.m),
    s.cons.hydro.total_energy_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.pressure(*s.m),
    s.prim.specific_internal_energy(*s.m),
    s.prim.sound_speed(*s.m),
    s.eos);
  auto lmax_f = sc.execute<tasks::hydro::update_max_characteristic_speed<D>>(
    flecsi::exec::on,
    *s.m,
    s.cons.hydro.mass_density(*s.m),
    s.prim.velocity(*s.m),
    s.prim.sound_speed(*s.m));
  s.dtmin_ = reduce<tasks::hydro::update_dtmin<D>, exec::fold::min>(
    flecsi::exec::on, *s.m, lmax_f);

  sc.execute<tasks::apply_boundaries<D>>(flecsi::exec::on,
    *s.m,
    s.icst.bmap(*s.gt),
    std::vector{s.cons.hydro.mass_density(*s.m),
      s.prim.pressure(*s.m),
      s.prim.specific_internal_energy(*s.m),
      s.cons.hydro.total_energy_density(*s.m)},
    std::vector{s.prim.velocity(*s.m), s.cons.hydro.momentum_density(*s.m)});

  /*--------------------------------------------------------------------------*
    Initialize time to 0
   *--------------------------------------------------------------------------*/
  sc.execute<tasks::init::init_time>(
    flecsi::exec::on, s.t(*s.gt), config["t0"].as<double>());

} // initialize

inline control<state, 1>::action<initialize<1>, cp::initialize> init_1d;
inline control<state, 2>::action<initialize<2>, cp::initialize> init_2d;
inline control<state, 3>::action<initialize<3>, cp::initialize> init_3d;

} // namespace hard
