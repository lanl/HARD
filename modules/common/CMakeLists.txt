# modules/common/CMakeLists.txt
#------------------------------------------------------------------------------#
# Common module (header-only library).
#   Provides:
#     library target "libcommon"
#   Headers (explicit list, no globbing).
#------------------------------------------------------------------------------#

add_library(libcommon INTERFACE)
add_library(libcommon::libcommon ALIAS libcommon)

set(modules_headers
  utils.hh
  tasks/utils.hh
)

target_sources(libcommon
  INTERFACE
    ${modules_headers}
)

target_include_directories(libcommon
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Keep these only if common headers use them in inline code.
target_link_libraries(libcommon
  INTERFACE
    yaml-cpp::yaml-cpp
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(libcommon
  INTERFACE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_CATALYST}>:ENABLE_CATALYST>
)

if(ENABLE_CATALYST)
  target_link_libraries(libcommon INTERFACE catalyst::catalyst)
endif()
