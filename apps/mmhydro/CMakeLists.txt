# hydro/CMakeLists.txt
#------------------------------------------------------------------------------#
# MMHydro application (executable).
#------------------------------------------------------------------------------#

add_executable(mmhydro
  mmhydro.cc
  actions/init.cc
  actions/rhs.cc
  actions/rk_stage_1.cc
  actions/rk_stage_2.cc
  actions/update_ls.cc
  actions/update_dt.cc
  actions/finalize.cc
  actions/analyze.cc
)

target_include_directories(mmhydro
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(mmhydro
  PRIVATE
    libcommon::libcommon
    libhydro::libhydro
    libgfm::libgfm
    libspec::libspec
    pybind11::pybind11
    singularity-eos::singularity-eos
    FleCSI::FleCSI
)

target_compile_definitions(mmhydro
  PRIVATE
    $<$<BOOL:${HARD_BENCHMARK_MODE}>:HARD_BENCHMARK_MODE>
    $<$<BOOL:${HARD_ENABLE_LEGION_TRACING}>:HARD_ENABLE_LEGION_TRACING>
    $<$<BOOL:${HARD_WRITE_CONTROL_INFO}>:HARD_WRITE_CONTROL_INFO>
)

if(ENABLE_CATALYST)
  target_link_libraries(hydro PRIVATE catalyst::catalyst)
  target_compile_definitions(hydro PRIVATE ENABLE_CATALYST)
endif()

if(ENABLE_UNIT_TESTS)

  set(FLECSI_UNIT_MAIN ${CMAKE_SOURCE_DIR}/apps/mmhydro/tests/main.cc)

endif()
install(TARGETS mmhydro)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/configs/
  DESTINATION share/HARD/mmhydro/configs
  FILES_MATCHING PATTERN "*.py"
  PATTERN "ci_configs" EXCLUDE
)
