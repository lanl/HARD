# modules/hydro/CMakeLists.txt
#------------------------------------------------------------------------------#
# Hydro module (header-only library).
#   Provides:
#     library target "libhydro"
#   Headers (explicit list, no globbing).
#------------------------------------------------------------------------------#

add_library(libhydro INTERFACE)
add_library(libhydro::libhydro ALIAS libhydro)

# You asked for an explicit list pattern; keep it as a set() so it is deterministic.
set(modules_headers
  constants.hh
  state.hh
  types.hh
  numerical_algorithms/riemann_solvers.hh
  numerical_algorithms/time_stepper.hh
  tasks/cons2prim.hh
  tasks/external_source.hh
  tasks/init.hh
  tasks/rhs.hh
  tasks/interface_fluxes.hh
  tasks/reconstruct.hh
  tasks/maxcharspeed.hh
  tasks/utils.hh
  tasks/initial_data/acoustic_wave.hh
  tasks/initial_data/kelvin_helm.hh
  tasks/initial_data/lw_implosion.hh
  tasks/initial_data/rayleigh_taylor.hh
  tasks/initial_data/sedov.hh
  tasks/initial_data/shock_tube.hh
)

target_sources(libhydro
  INTERFACE
    ${modules_headers}
)

target_include_directories(libhydro
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# libhydro uses the shared spec module + external deps.
target_link_libraries(libhydro
  INTERFACE
    libcommon::libcommon
    libspec::libspec
    yaml-cpp::yaml-cpp
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(libhydro
  INTERFACE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_RADIATION}>:ENABLE_RADIATION>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
    $<$<BOOL:${ENABLE_FLECSOLVE}>:ENABLE_FLECSOLVE>
)

if(ENABLE_CATALYST)
  target_link_libraries(libhydro INTERFACE catalyst::catalyst)
endif()

if(ENABLE_FLECSOLVE)
  target_link_libraries(libhydro INTERFACE flecsolve::flecsolve)
endif()
