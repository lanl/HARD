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
  io.hh
  labels.hh
  limiter.hh
  mesh.hh
  runtime.hh
  types.hh
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
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/modules>
    $<INSTALL_INTERFACE:include>   # if you install later
)

# Spec headers may depend on these (keep only what headers actually use).
target_link_libraries(libspec
  INTERFACE
    singularity-eos::singularity-eos
    FleCSI::FleCSI
    pybind11::headers
    Python::Python
)

target_compile_definitions(libspec
  INTERFACE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
)

if(ENABLE_CATALYST)
  target_link_libraries(libspec INTERFACE catalyst::catalyst)
endif()
