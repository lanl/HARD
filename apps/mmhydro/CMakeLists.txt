# mmhydro/CMakeLists.txt
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
    libspec::libspec
    singularity-eos::singularity-eos
    FleCSI::FleCSI
    pybind11::pybind11
    Python::Python
    libgfm::libgfm
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

if(ENABLE_VERIFICATION AND ENABLE_UNIT_TESTS AND NOT ENABLE_CATALYST)

  set(TEST_VERIFY_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_verify.sh")
  set(TEST_DISTRIBUTED_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_distributed.sh")
  set(TEST_CONVERGENCE_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_convergence.sh")

  # Accuracy tests
  add_test(NAME mmhydro_accuracy_sod_1d
    COMMAND ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/mmhydro/mmhydro
    ${CMAKE_SOURCE_DIR}/configs/sod_mm.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

endif()

install(TARGETS mmhydro)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/configs/
  DESTINATION share/HARD/mmhydro/configs
  FILES_MATCHING PATTERN "*.py"
  PATTERN "ci_configs" EXCLUDE
)
