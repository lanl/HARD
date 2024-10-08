#include "advance.hh"
#include "analyze.hh"
#include "finalize.hh"
#include "init.hh"
#include "options.hh"
#include "state.hh"

#include <spec/runtime.hh>

#include <flecsi/runtime.hh>
#include <yaml-cpp/yaml.h>

#ifdef USE_CATALYST
#include "catalyst/adaptor.hh"
#include "tasks/catalyst.hh"
#endif

#include <flecsi/util/unit.hh>

using namespace flecsi;
using namespace spec;
using namespace hard;

int
main(int argc, char ** argv) {
  // Output the control model and actions if the user has enabled
  // HARD_WRITE_CONTROL_INFO.
#if defined(HARD_WRITE_CONTROL_INFO)
  control<state, 3>::write_graph("HARD", "cm.dot");
  control<state, 3>::write_actions("HARD", "actions.dot");
#endif

  const flecsi::getopt g;
  try {
    g(argc, argv);
  }
  catch(const std::logic_error & e) {
    std::cerr << e.what() << '\n' << g.usage(argc ? argv[0] : "");
    return 1;
  }

  const flecsi::run::dependencies_guard dg;
  run::config cfg;
  cfg.flog.tags = {opt::flog_tags};
  cfg.flog.verbose = {opt::flog_verbose};
  cfg.flog.process = {opt::flog_process};

#if FLECSI_BACKEND == FLECSI_BACKEND_legion
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  cfg.legion = {"", "-ll:gpu", "1"};
#elif defined(KOKKOS_ENABLE_OPENMP) && defined(REALM_USE_OPENMP)
  cfg.legion = {"", "-ll:ocpu", "1", "-ll:onuma", "0"};
#endif
#endif

  runtime run{cfg};
  flog::add_output_stream("clog", std::clog, true);

  YAML::Node config = YAML::LoadFile(opt::config.value());

#ifdef HARD_BENCHMARK_MODE
  // Write header if process 0 and header==1
  if(flecsi::process() == 0 && opt::header.value()) {
    std::ofstream runtime_file;
    runtime_file.open("result/runtimes.txt", std::ios_base::app);
    runtime_file << "process;threads;max_iterations;total_runtime;iteration_"
                    "runtimes;dimension="
                 << opt::dimension.value() << ";\n";
    runtime_file.close();
  }
#endif

  return dispatch<control, state>(run,
    opt::dimension.value(),
    config["t0"].as<double>(),
    config["tf"].as<double>(),
    config["max_steps"].as<std::size_t>(),
    config["cfl"].as<double>(),
    config["max_dt"].as<double>(),
    config["log_frequency"].as<std::size_t>(),
    config["output_frequency"].as<std::size_t>());
} // main
