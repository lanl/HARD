# hydro/CMakeLists.txt
#------------------------------------------------------------------------------#
# Hydro application (executable).
#------------------------------------------------------------------------------#

add_executable(hydro
  hydro.cc
)

target_include_directories(hydro
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(hydro
  PRIVATE
    libcommon::libcommon
    libhydro::libhydro
    libspec::libspec
    yaml-cpp::yaml-cpp
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(hydro
  PRIVATE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_RADIATION}>:ENABLE_RADIATION>
)

if(ENABLE_CATALYST)
  target_link_libraries(hydro PRIVATE catalyst::catalyst)
  target_compile_definitions(hydro PRIVATE ENABLE_CATALYST)
endif()

if(ENABLE_FLECSOLVE)
  target_link_libraries(hydro PRIVATE flecsolve::flecsolve)
  target_compile_definitions(hydro PRIVATE ENABLE_FLECSOLVE)
endif()
