#ifndef HARD_HYDRO_OPTIONS_HH
#define HARD_HYDRO_OPTIONS_HH

#include <flecsi/execution.hh>

namespace hard::opt {

inline flecsi::program_option<std::string> config("python file",
  "The python config file.",
  1,
  [](std::string const & v, std::stringstream & ss) {
    return v.find(".py") != std::string::npos
             ? true
             : (ss << "file(" << v << ") has invalid suffix") && false;
  });

inline flecsi::program_option<unsigned int> dimension("Hard Options",
  "dimension,d",
  "Specify the dimension of the solver (default: 3).",
  {{flecsi::option_default, 3}});

inline flecsi::program_option<unsigned int> resolution("Hard Options",
  "resolution,r",
  "Specify the lowest level resolution for the grid.",
  {{flecsi::option_default, 0}});

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

#endif // HARD_HYDRO_OPTIONS_HH
