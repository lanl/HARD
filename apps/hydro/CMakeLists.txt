# hydro/CMakeLists.txt
#------------------------------------------------------------------------------#
# Hydro application (executable).
#------------------------------------------------------------------------------#

add_executable(hydro
  hydro.cc
  actions/init.cc
  actions/rhs.cc
  actions/rk_stage_1.cc
  actions/rk_stage_2.cc
  actions/update_dt.cc
  actions/finalize.cc
  actions/analyze.cc
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
    singularity-eos::singularity-eos
    FleCSI::FleCSI
    pybind11::pybind11
    Python::Python
)

target_compile_definitions(hydro
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

  set(FLECSI_UNIT_MAIN ${CMAKE_SOURCE_DIR}/apps/hydro/tests/main.cc)

  flecsi_add_test(mesh
    SOURCES
    tests/mesh.cc
  )

endif()

if(ENABLE_VERIFICATION AND ENABLE_UNIT_TESTS AND NOT ENABLE_CATALYST)

  set(TEST_VERIFY_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_verify.sh")
  set(TEST_DISTRIBUTED_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_distributed.sh")
  set(TEST_CONVERGENCE_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_convergence.sh")

  file(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/tests/rk_3d
    ${CMAKE_BINARY_DIR}/tests/sod_3d
  )

  # Accuracy tests
  add_test(NAME accuracy_rk_1d
    COMMAND ${TEST_VERIFY_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/rk.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_sod_1d
    COMMAND ${TEST_VERIFY_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/sod.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_leblanc_1d
    COMMAND ${TEST_VERIFY_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/leblanc.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  # add_test(NAME accuracy_sedov_1d
  #   COMMAND ${TEST_VERIFY_SCRIPT}
  #   ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
  #   ${CMAKE_SOURCE_DIR}/configs/ci_configs/sedov.py -d1
  #   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  # )
  add_test(NAME accuracy_acoustic_1d
    COMMAND ${TEST_VERIFY_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_acoustic_2d
    COMMAND ${TEST_VERIFY_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave.py -d2
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

  # Error with options for HPX for ranks > 1
  if(NOT FleCSI_BACKEND STREQUAL "hpx")
    # Distributed tests
    add_test(NAME hard-distributed
      COMMAND ${TEST_DISTRIBUTED_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/distr-rk.py -d1)
  endif()

  # Convergence tests
  add_test(NAME total_convergence_testing
    COMMAND ${TEST_CONVERGENCE_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    6 10
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave_total_convergence.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME spatial_convergence_testing
    COMMAND ${TEST_CONVERGENCE_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/hydro/hydro
    4 8
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave_spatial_convergence.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

endif()

install(TARGETS hydro)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/configs/
  DESTINATION share/HARD/hydro/configs
  FILES_MATCHING PATTERN "*.py"
  PATTERN "ci_configs" EXCLUDE
)
