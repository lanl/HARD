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

if(ENABLE_UNIT_TESTS)

  set(FLECSI_UNIT_MAIN ${CMAKE_SOURCE_DIR}/hydro/tests/main.cc)

  flecsi_add_test(mesh
    SOURCES
    tests/mesh.cc
  )

endif()

if(ENABLE_VERIFICATION AND ENABLE_UNIT_TESTS)

  if(ENABLE_UNIT_TESTS AND NOT ENABLE_CATALYST)

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
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/rk.yaml -d1
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    add_test(NAME accuracy_sod_1d
      COMMAND ${TEST_VERIFY_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/sod.yaml -d1
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    add_test(NAME accuracy_leblanc_1d
      COMMAND ${TEST_VERIFY_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/leblanc.yaml -d1
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    # add_test(NAME accuracy_sedov_1d
    #   COMMAND ${TEST_VERIFY_SCRIPT}
    #   ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
    #   ${CMAKE_SOURCE_DIR}/configs/ci_configs/sedov.yaml -d1
    #   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # )
    add_test(NAME accuracy_acoustic_1d
      COMMAND ${TEST_VERIFY_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave.yaml -d1
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    add_test(NAME accuracy_acoustic_2d
      COMMAND ${TEST_VERIFY_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave.yaml -d2
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    # Distributed tests
    add_test(NAME hard-distributed
      COMMAND ${TEST_DISTRIBUTED_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/distr-rk.yaml -d1)

    # Convergence tests
    add_test(NAME convergence_testing
      COMMAND ${TEST_CONVERGENCE_SCRIPT}
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/hydro/hydro
      6 10
      ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave_convergence.yaml -d1
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  endif()
endif()

install(TARGETS hydro)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/configs/
  DESTINATION share/HARD/hydro/configs
  FILES_MATCHING PATTERN "*.yaml"
  PATTERN "ci_configs" EXCLUDE
)
