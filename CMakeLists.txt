#------------------------------------------------------------------------------#
# Copyright (c) 2023, Triad National Security, LLC.
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.19)

#------------------------------------------------------------------------------#
# Add local module path.
#------------------------------------------------------------------------------#

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------#
# CMake policies.
#------------------------------------------------------------------------------#

cmake_policy(SET CMP0144 NEW) # ``PACKAGE_NAME''_ROOT

#------------------------------------------------------------------------------#
# Project.
#------------------------------------------------------------------------------#

project(HARD LANGUAGES C CXX)
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS
    "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------#
# Dependencies (load early).
#------------------------------------------------------------------------------#

find_package(FleCSI 2 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(singularity-eos REQUIRED COMPONENTS Interface)
find_package(Eigen3 REQUIRED)
find_package(flecsolve REQUIRED)


#------------------------------------------------------------------------------#
# Catalyst
#------------------------------------------------------------------------------#

option(ENABLE_CATALYST "Enable Catalyst output" OFF)
if(ENABLE_CATALYST)
  find_package(catalyst REQUIRED)
endif()

#------------------------------------------------------------------------------#
# Control-Model Output
#------------------------------------------------------------------------------#

option(HARD_WRITE_CONTROL_INFO
  "Output the control model graph and actions at startup" OFF)
mark_as_advanced(HARD_WRITE_CONTROL_INFO)

#------------------------------------------------------------------------------#
# Benchmark Mode
#------------------------------------------------------------------------------#

option(HARD_BENCHMARK_MODE
  "Benchmark mode disables I/O and adds time measurement" OFF)
mark_as_advanced(HARD_BENCHMARK_MODE)

option(HARD_ENABLE_LEGION_TRACING
  "Enable Legion Tracing" OFF)
mark_as_advanced(HARD_ENABLE_LEGION_TRACING)

#------------------------------------------------------------------------------#
# Formatting
#------------------------------------------------------------------------------#

option(ENABLE_FORMAT "Enable format target" OFF)
mark_as_advanced(ENABLE_FORMAT)

if(ENABLE_FORMAT)
  include(format)

  set(CLANG_FORMAT_VERSION "20...<21" CACHE STRING
      "Set the required version (major[.minor[.patch]]) of clang-format")
  mark_as_advanced(CLANG_FORMAT_VERSION)

  flecsi_add_format_target(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}
    "${CLANG_FORMAT_VERSION}")
endif()

#------------------------------------------------------------------------------#
# Documentation (kept as-is; optional).
#------------------------------------------------------------------------------#

option(ENABLE_DOCUMENTATION "Enable documentation" OFF)
mark_as_advanced(ENABLE_DOCUMENTATION)

if(ENABLE_DOCUMENTATION)
  find_package(LATEX REQUIRED)
  find_program(PDFTOPPM_EXECUTABLE pdftoppm PATHS /usr/bin)
  add_custom_target(structure
    COMMAND
      ${PDFLATEX_COMPILER} ${PROJECT_SOURCE_DIR}/doc/tex/structure.tex &&
      ${PDFTOPPM_EXECUTABLE} -jpeg structure.pdf > structure.jpg
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

  include(FleCSI/sphinx)

  file(COPY ${CMAKE_SOURCE_DIR}/doc/sphinx
    DESTINATION ${CMAKE_BINARY_DIR}
  )
  file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/.sphinx)
  file(RENAME ${CMAKE_BINARY_DIR}/sphinx
    ${CMAKE_BINARY_DIR}/.sphinx
  )
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/.sphinx/hard/src)
  execute_process(COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/spec
    ${CMAKE_BINARY_DIR}/.sphinx/hard/src/spec
  )
  execute_process(COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/app
    ${CMAKE_BINARY_DIR}/.sphinx/hard/src/app
  )

  flecsi_add_sphinx_target(main
    CONFIG ${CMAKE_BINARY_DIR}/.sphinx
    OUTPUT ${CMAKE_BINARY_DIR}/doc
  )

  find_package(Git REQUIRED)
  add_custom_target(deploy-doc
    COMMAND
      make doc &&
      echo "Updating pages" &&
      ([ -e pages ] ||
        ${GIT_EXECUTABLE} clone -q --single-branch --branch pages
        git@gitlab.lanl.gov:cdss24/hard.git pages &&
        cd pages &&
        ${GIT_EXECUTABLE} rm -qr . && ${GIT_EXECUTABLE} reset -q &&
        ${GIT_EXECUTABLE} checkout .gitignore &&
        ${GIT_EXECUTABLE} checkout .gitlab-ci.yml) &&
      echo "Updating pages" &&
        cp -rT doc pages/public &&
        (cd pages && ${GIT_EXECUTABLE} add -A .) &&
      echo "Updated pages are in ${CMAKE_BINARY_DIR}/pages"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

#------------------------------------------------------------------------------#
# Testing 
#------------------------------------------------------------------------------#

option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)
mark_as_advanced(ENABLE_UNIT_TESTS)

include(FleCSI/unit)

if(ENABLE_UNIT_TESTS)
  flecsi_enable_testing()
endif()

#------------------------------------------------------------------------------#
# modules and apps.
#------------------------------------------------------------------------------#

add_subdirectory(modules/common)
add_subdirectory(modules/spec)
add_subdirectory(modules/hydro)
add_subdirectory(modules/rad)
add_subdirectory(apps/hydro)
add_subdirectory(apps/radhydro)

#------------------------------------------------------------------------------#
# Optional verification
#------------------------------------------------------------------------------#

option(ENABLE_VERIFICATION "Enable exactpack verification" OFF)
mark_as_advanced(ENABLE_FORMAT)

if(ENABLE_VERIFICATION AND ENABLE_UNIT_TESTS)

  find_package(Python3 REQUIRED COMPONENTS Interpreter)

  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import exactpack"
    RESULT_VARIABLE EXACTPACK_FOUND
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(NOT EXACTPACK_FOUND EQUAL 0)
    message(FATAL_ERROR "ExactPack Python module not found. Please install it (e.g., with 'pip install exactpack').")
  endif()

  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import matplotlib"
    RESULT_VARIABLE MATPLOTLIB_FOUND
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(NOT MATPLOTLIB_FOUND EQUAL 0)
    message(FATAL_ERROR "Matplotlib Python module not found. Please install it (e.g., with 'pip install matplotlib').")
  endif()

  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import yaml"
    RESULT_VARIABLE PYYAML_FOUND
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(NOT PYYAML_FOUND EQUAL 0)
    message(FATAL_ERROR "Yaml Python module not found. Please install it (e.g., with 'pip install pyyaml').")
  endif()

  set(CHECK_SCRIPT ${CMAKE_SOURCE_DIR}/tools/exactpack_verify.py)
  set(CONVERGENCE_SCRIPT ${CMAKE_SOURCE_DIR}/tools/convergence_verify.py)

  if(CONFIG STREQUAL "")
    set(CONFIG sod.yaml)
  endif()

  add_custom_target(verify
    COMMAND ${Python3_EXECUTABLE} ${CHECK_SCRIPT} \$\(CONFIG\) \$\(DIMENSION\) "--plot"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
  )

  add_custom_target(convergence
    COMMAND ${Python3_EXECUTABLE} ${CONVERGENCE_SCRIPT} \$\(CONFIG\) \$\(DIMENSION\)
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
  )
endif()
