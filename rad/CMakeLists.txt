# rad/CMakeLists.txt
#------------------------------------------------------------------------------#
# Rad application (executable).
#------------------------------------------------------------------------------#

add_executable(rad
  rad.cc
  actions/analyze.cc
  actions/init.cc
  actions/finalize.cc
  actions/rk_stage_1.cc
  actions/rk_stage_2.cc
  actions/couple_hydro_radiation.cc
  actions/radiation.cc
  actions/update_dt.cc
  actions/rhs.cc
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

if(ENABLE_CATALYST)
  target_link_libraries(hydro PRIVATE catalyst::catalyst)
  target_compile_definitions(hydro PRIVATE ENABLE_CATALYST)
endif()

install(TARGETS rad)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/configs/
  DESTINATION share/HARD/rad/configs
  FILES_MATCHING PATTERN "*.yaml"
  PATTERN "ci_configs" EXCLUDE
)
