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
  tasks/init.hh
  tasks/rad.hh
  tasks/rad_root.hh
  state.hh
  utils.hh
  constants.hh
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
    $<$<BOOL:${ENABLE_RADIATION}>:ENABLE_RADIATION>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
    $<$<BOOL:${ENABLE_FLECSOLVE}>:ENABLE_FLECSOLVE>
)

target_link_libraries(librad INTERFACE flecsolve::flecsolve)
