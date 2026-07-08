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


set(MMHYDRO_HPX_TEST_ENV
  "HPX_SMALL_STACK_SIZE=0x400000"
  "HPX_MEDIUM_STACK_SIZE=0x800000"
  "HPX_LARGE_STACK_SIZE=0x1000000"
)

if(ENABLE_UNIT_TESTS AND (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP))

  flecsi_add_test(mmhydro_level_set_test
    SOURCES
    tests/level_set.cc
    actions/update_ls.cc
  )

  target_include_directories(mmhydro_level_set_test
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests
      ${CMAKE_SOURCE_DIR}/modules
  )

  target_link_libraries(mmhydro_level_set_test
    PRIVATE
      libcommon::libcommon
      libspec::libspec
      libgfm::libgfm
      FleCSI::FleCSI
  )

  target_compile_definitions(mmhydro_level_set_test
    PRIVATE
      IN_TEST
  )

endif()

if(ENABLE_VERIFICATION AND ENABLE_UNIT_TESTS)

  set(TEST_VERIFY_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_verify.sh")
  set(TEST_DISTRIBUTED_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_distributed.sh")
  set(TEST_CONVERGENCE_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_convergence.sh")
  set(TEST_CONVERGENCE_MM_SCRIPT "${CMAKE_SOURCE_DIR}/tools/test_convergence_mm.sh")

  set(LEGION_GPU_FLAGS "")
  if(ENABLE_LEGION AND ENABLE_CUDA)
    set(LEGION_GPU_FLAGS "-ll:gpu 1 -ll:show_rsrv")
  endif()

  # Accuracy tests
  add_test(NAME mmhydro_accuracy_sod_1d
    COMMAND ${MPIEXEC_EXECUTABLE} -n 1 ${CMAKE_BINARY_DIR}/apps/mmhydro/mmhydro
    ${LEGION_GPU_FLAGS} ${CMAKE_SOURCE_DIR}/configs/ci_configs/sod_mm.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  set_tests_properties(mmhydro_accuracy_sod_1d PROPERTIES
    ENVIRONMENT "${MMHYDRO_HPX_TEST_ENV}"
  )

  if(ENABLE_UNIT_TESTS AND (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP))

  add_test(NAME mmhydro_contact_move_convergence_1d
    COMMAND ${TEST_CONVERGENCE_MM_SCRIPT}
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/apps/mmhydro/mmhydro
    7 11
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/contact_move_convergence.py -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  set_tests_properties(mmhydro_contact_move_convergence_1d PROPERTIES
    ENVIRONMENT "${MMHYDRO_HPX_TEST_ENV}"
  )
  endif()

endif()

install(TARGETS mmhydro)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/configs/
  DESTINATION share/HARD/mmhydro/configs
  FILES_MATCHING PATTERN "*.py"
  PATTERN "ci_configs" EXCLUDE
)
