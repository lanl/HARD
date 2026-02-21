# modules/spec/CMakeLists.txt
#------------------------------------------------------------------------------#
# Spec module (header-only library).
#   Provides:
#     library target "libspec"
#   Headers (explicit list, no globbing).
#------------------------------------------------------------------------------#

add_library(libspec INTERFACE)
add_library(libspec::libspec ALIAS libspec)

# Explicit list (and example of how to keep it readable).
set(modules_headers
  control.hh
  eos.hh
  labels.hh
  limiter.hh
  mesh.hh
  rad_function.hh
  runtime.hh
  types.hh
  utils.hh
  tasks/boundaries/boundary.hh
  tasks/boundaries/dirichlet.hh
  tasks/boundaries/flow.hh
  tasks/boundaries/reflective.hh
  tasks/io.hh

)

# Add headers explicitly so IDEs show them and installs/exports can reference them.
target_sources(libspec
  INTERFACE
    ${modules_headers}
)

target_include_directories(libspec
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Spec headers may depend on these (keep only what headers actually use).
target_link_libraries(libspec
  INTERFACE
    yaml-cpp::yaml-cpp
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(libspec
  INTERFACE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_RADIATION}>:ENABLE_RADIATION>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
    $<$<BOOL:${ENABLE_FLECSOLVE}>:ENABLE_FLECSOLVE>
)

if(ENABLE_CATALYST)
  target_link_libraries(libspec INTERFACE catalyst::catalyst)
endif()

if(ENABLE_FLECSOLVE)
  target_link_libraries(libspec INTERFACE flecsolve::flecsolve)
endif()
