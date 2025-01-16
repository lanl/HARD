#ifndef HARD_OPTIONS_HH
#define HARD_OPTIONS_HH

#include <flecsi/execution.hh>

namespace hard::opt {

inline flecsi::program_option<std::string> config("yaml file",
  "The yaml config file.",
  1,
  [](std::string const & v, std::stringstream & ss) {
    return v.find(".yaml") != std::string::npos
             ? true
             : (ss << "file(" << v << ") has invalid suffix") && false;
  });

#ifdef USE_CATALYST
inline flecsi::program_option<std::string> catalyst_script("catalyst script",
  "The catalyst script file.",
  1,
  [](flecsi::any const & v, std::stringstream & ss) {
    const std::string value = flecsi::option_value<std::string>(v);
    return value.find(".py") != std::string::npos
             ? true
             : (ss << "file(" << value << ") has invalid suffix") && false;
  });
#endif

inline flecsi::program_option<unsigned int> dimension("Hard Options",
  "dimension,d",
  "Specify the dimension of the solver (default: 3).",
  {{flecsi::option_default, 3}});

inline flecsi::program_option<std::string> source_fds("Hard Options",
  "fds_file",
  "Specify where the file lives (default: source_fds.txt)",
  {{flecsi::option_default, "source_fds.txt"}});

inline flecsi::program_option<unsigned int> colors("Hard Options",
  "colors,c",
  "Specify the number of colors (default: num processes).",
  {{flecsi::option_default, 0}});

inline flecsi::program_option<std::string> flog_tags("FLOG Options",
  "tags,t",
  "Specify the flog tags to enable.",
  {{flecsi::option_default, "all"}});

inline flecsi::program_option<int> flog_verbose("FLOG Options",
  "verbose,v",
  "Enable verbose output. Passing '-1' will strip any additional"
  " decorations added by flog and will only output the user's message.",
  {{flecsi::option_default, 0}});

inline flecsi::program_option<int> flog_process("FLOG Options",
  "process,p",
  "Specify output process. Passing '-1' will enable output from all processes.",
  {{flecsi::option_default, 0}});
#ifdef HARD_BENCHMARK_MODE
inline flecsi::program_option<bool> header("Header",
  "header,h",
  "Write header line into result file. Passing '1' will write header line",
  {{flecsi::option_default, 0}});
#endif
} // namespace hard::opt

#endif // HARD_OPTIONS_HH
