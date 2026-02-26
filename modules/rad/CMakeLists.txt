# modules/rad/CMakeLists.txt
#------------------------------------------------------------------------------#
# Rad module (header-only library).
#   Provides:
#     library target "librad"
#------------------------------------------------------------------------------#

add_library(librad INTERFACE)
add_library(librad::librad ALIAS librad)

# You asked for an explicit list pattern; keep it as a set() so it is deterministic.
set(modules_headers
  state.hh
  constants.hh
  flecsolvers.hh
  types.hh
  linsolve.hh
  tasks/init.hh
  tasks/rad.hh
  tasks/rad_root.hh
  tasks/interface_fluxes.hh
  tasks/operators.hh
  tasks/time_derivative.hh
  tasks/utils.hh
  tasks/initial_data/all_initial_data.hh
  tasks/initial_data/heating_and_cooling_afld.hh
  tasks/initial_data/heating_and_cooling.hh
  tasks/initial_data/implosion_forced_T.hh
  tasks/initial_data/kelvin_helm_rad.hh
  tasks/initial_data/rad_rank_hugoniot.hh
  tasks/initial_data/richtmyer_meshkov.hh
)

target_sources(librad
  INTERFACE
    ${modules_headers}
)

target_include_directories(librad
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# libhydro uses the shared spec module + external deps.
target_link_libraries(librad
  INTERFACE
    libcommon::libcommon
    libspec::libspec
    yaml-cpp::yaml-cpp
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(librad
  INTERFACE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
)

target_link_libraries(librad INTERFACE flecsolve::flecsolve)
