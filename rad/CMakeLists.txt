# rad/CMakeLists.txt
#------------------------------------------------------------------------------#
# Rad application (executable).
#------------------------------------------------------------------------------#

add_executable(rad
  rad.cc
)

target_include_directories(rad
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(rad
  PRIVATE
    libcommon::libcommon
    libhydro::libhydro
    librad::librad
    libspec::libspec
    yaml-cpp::yaml-cpp
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(rad
  PRIVATE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
    $<$<BOOL:${ENABLE_RADIATION}>:ENABLE_RADIATION>
)
