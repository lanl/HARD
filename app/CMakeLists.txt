set(HARD_SOURCES
  hard.cc
  numerical_algorithms/root_finder.cc
)

if(ENABLE_CATALYST)
  set(HARD_SOURCES
    ${HARD_SOURCES}
    catalyst/types.cc
    catalyst/adaptor.cc
  ) 
endif()

add_executable(hard ${HARD_SOURCES})

if(Kokkos_ENABLE_CUDA AND CMAKE_BUILD_TYPE STREQUAL "Debug" AND
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 13.0.1)
  message(WARNING "Disabling GPU debuginfo for Clang > 13.0.1 due to "
                  "https://github.com/llvm/llvm-project/issues/58491")
  target_compile_options(hard PUBLIC -Xarch_device -g0)
endif()

target_link_libraries(hard
	Spec::Spec
	FleCSI::FleCSI
	yaml-cpp::yaml-cpp
  singularity-eos::singularity-eos
)


if(ENABLE_CATALYST)
  target_compile_definitions(hard
    PRIVATE
      USE_CATALYST=1)

  target_link_libraries(hard
    catalyst::catalyst)
endif()

if(HARD_WRITE_CONTROL_INFO)
  target_compile_definitions(hard PUBLIC HARD_WRITE_CONTROL_INFO)
endif()

if(HARD_BENCHMARK_MODE)
  target_compile_definitions(hard PUBLIC HARD_BENCHMARK_MODE)
endif()

if(ENABLE_RADIATION)
  target_compile_definitions(hard PUBLIC ENABLE_RADIATION)
endif()

if(HARD_ENABLE_LEGION_TRACING)
    target_compile_definitions(hard PUBLIC HARD_ENABLE_LEGION_TRACING)
endif()

if(ENABLE_UNIT_TESTS AND NOT ENABLE_CATALYST)
  file(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/tests/rk_3d
    ${CMAKE_BINARY_DIR}/tests/sod_3d
  )

if(NOT ENABLE_RADIATION)
  # Accuracy tests 
  add_test(NAME accuracy_rk_1d
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_verify.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/rk.yaml -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_sod_1d
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_verify.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/sod.yaml -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_leblanc_1d
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_verify.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/leblanc.yaml -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_sedov_1d
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_verify.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/sedov.yaml -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_test(NAME accuracy_acoustic_1d
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_verify.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave.yaml -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

if(NOT ENABLE_RADIATION)
  # Distributed tests 
  add_test(NAME hard-distributed
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_distributed.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/distr-rk.yaml -d1)
endif()

if(NOT ENABLE_RADIATION)
  # Convergence tests
  add_test(NAME convergence_testing
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_convergence.sh
    ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/hard
    7 12
    ${CMAKE_SOURCE_DIR}/configs/ci_configs/acoustic-wave_convergence.yaml -d1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

endif()

install(TARGETS hard)
