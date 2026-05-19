# modules/gfm/CMakeLists.txt
#------------------------------------------------------------------------------#
# GFM module (header-only library).
#   Provides:
#     library target "libgfm"
#   Headers (explicit list, no globbing).
#------------------------------------------------------------------------------#

add_library(libgfm INTERFACE)
add_library(libgfm::libgfm ALIAS libgfm)

# You asked for an explicit list pattern; keep it as a set() so it is deterministic.
set(modules_headers
  types.hh
  state.hh
  tasks/ls_core.hh
  tasks/extrapolation.hh
  tasks/reconstruct.hh
  tasks/interface_fluxes_gfm.hh
  tasks/cons2prim.hh
  tasks/initial_data/shock_tube_mm.hh
)

target_sources(libgfm
  INTERFACE
    ${modules_headers}
)

target_include_directories(libgfm
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/modules>
    $<INSTALL_INTERFACE:include>   # if you install later
)

# libgfm uses the shared spec module + external deps.
target_link_libraries(libgfm
  INTERFACE
    libcommon::libcommon
    libspec::libspec
    pybind11::headers
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(libgfm
  INTERFACE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
)

if(ENABLE_CATALYST)
  target_link_libraries(libgfm INTERFACE catalyst::catalyst)
endif()
