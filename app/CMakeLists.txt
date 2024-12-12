set(FLASTRO_SOURCES
  flastro.cc
  numerical_algorithms/root_finder.cc
)

if(ENABLE_CATALYST)
  set(FLASTRO_SOURCES
    ${FLASTRO_SOURCES}
    catalyst/types.cc
    catalyst/adaptor.cc
  ) 
endif()

add_executable(flastro ${FLASTRO_SOURCES})

if(Kokkos_ENABLE_CUDA AND CMAKE_BUILD_TYPE STREQUAL "Debug" AND
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 13.0.1)
  message(WARNING "Disabling GPU debuginfo for Clang > 13.0.1 due to "
                  "https://github.com/llvm/llvm-project/issues/58491")
  target_compile_options(flastro PUBLIC -Xarch_device -g0)
endif()

target_link_libraries(flastro
	Spec::Spec
	FleCSI::FleCSI
	yaml-cpp
)


if(ENABLE_CATALYST)
  target_compile_definitions(flastro
    PRIVATE
      USE_CATALYST=1)

  target_link_libraries(flastro
    catalyst::catalyst)
endif()

if(FLASTRO_WRITE_CONTROL_INFO)
  target_compile_definitions(flastro PUBLIC FLASTRO_WRITE_CONTROL_INFO)
endif()

if(FLASTRO_BENCHMARK_MODE)
  target_compile_definitions(flastro PUBLIC FLASTRO_BENCHMARK_MODE)
endif()

if(DISABLE_RADIATION)
  target_compile_definitions(flastro PUBLIC DISABLE_RADIATION)
endif()

if(FLASTRO_ENABLE_LEGION_TRACING)
    target_compile_definitions(flastro PUBLIC FLASTRO_ENABLE_LEGION_TRACING)
endif()

if(ENABLE_UNIT_TESTS AND NOT ENABLE_CATALYST)
  file(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/tests/rk_3d
    ${CMAKE_BINARY_DIR}/tests/sod_3d
  )

  add_test(NAME flastro_rk_3d
    COMMAND flastro ${CMAKE_SOURCE_DIR}/configs/rk.yaml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/rk_3d)
  add_test(NAME flastro_sod_3d
    COMMAND flastro ${CMAKE_SOURCE_DIR}/configs/sod.yaml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/sod_3d)

  add_test(NAME flastro-distributed
    COMMAND ${CMAKE_SOURCE_DIR}/app/tests/test_distributed.sh 
      ${MPIEXEC_EXECUTABLE} ${CMAKE_BINARY_DIR}/app/flastro 
      ${CMAKE_SOURCE_DIR}/configs/distr-rk.yaml -d1)

endif()
