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

set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS
    "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------#
# Find the FleCSI package. This must be done early so that the FleCSI cmake
# packages will be available to the project, e.g., 'include(FleCSI/sphinx)'.
#
# CMake documentation for FleCSI client applications is included in the
# FleCSI userguide.
#------------------------------------------------------------------------------#

if(NOT FleCSI_FOUND)
  find_package(FleCSI 2 REQUIRED)
endif()

#------------------------------------------------------------------------------#
# Yaml
#------------------------------------------------------------------------------#

find_package(yaml-cpp REQUIRED)

#------------------------------------------------------------------------------#
# Singularity
#------------------------------------------------------------------------------#

find_package(singularity-eos REQUIRED)

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

#------------------------------------------------------------------------------#
# Benchmark Mode
#------------------------------------------------------------------------------#

option(HARD_ENABLE_LEGION_TRACING
  "Enable Legion Tracing" OFF)
mark_as_advanced(HARD_ENABLE_LEGION_TRACING)

#------------------------------------------------------------------------------#
# Disable Radiation
#------------------------------------------------------------------------------#

option(DISABLE_RADIATION
  "Disables radiation in evolution" OFF)
mark_as_advanced(DISABLE_RADIATION)

#------------------------------------------------------------------------------#
# Formatting
#------------------------------------------------------------------------------#

option(ENABLE_FORMAT "Enable format target" OFF)
mark_as_advanced(ENABLE_FORMAT)

if(ENABLE_FORMAT)
  include(format)

  set(CLANG_FORMAT_VERSION "13...<14" CACHE STRING
      "Set the required version (major[.minor[.patch]]) of clang-format")
  mark_as_advanced(CLANG_FORMAT_VERSION)

  flecsi_add_format_target(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}
    "${CLANG_FORMAT_VERSION}")
endif()

#------------------------------------------------------------------------------#
# Documentation
#------------------------------------------------------------------------------#

option(ENABLE_DOCUMENTATION "Enable documentation" OFF)
mark_as_advanced(ENABLE_DOCUMENTATION)

if(ENABLE_DOCUMENTATION)
  #----------------------------------------------------------------------------#
  # README.md figure
  #----------------------------------------------------------------------------#

  find_package(LATEX REQUIRED)
  find_program(PDFTOPPM_EXECUTABLE pdftoppm PATHS /usr/bin)
  add_custom_target(structure
    COMMAND
      ${PDFLATEX_COMPILER} ${PROJECT_SOURCE_DIR}/doc/tex/structure.tex &&
      ${PDFTOPPM_EXECUTABLE} -jpeg structure.pdf > structure.jpg
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

  #----------------------------------------------------------------------------#
  # Sphinx & Doxygen
  #----------------------------------------------------------------------------#

  include(FleCSI/doxygen)
  include(FleCSI/sphinx)

  # Avoid race condition by proactively creating output directory
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/.doxygen/xml)
  flecsi_add_doxygen_target(hard
    CONFIGS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/hard.conf.in
      ${CMAKE_SOURCE_DIR}/doc/doxygen/conf.in
  )

  file(COPY ${CMAKE_SOURCE_DIR}/doc/sphinx
    DESTINATION ${CMAKE_BINARY_DIR}
  )
  file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/.sphinx)
  file(RENAME ${CMAKE_BINARY_DIR}/sphinx
    ${CMAKE_BINARY_DIR}/.sphinx
  )
  # Create output directory
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

if(ENABLE_UNIT_TESTS)
  include(FleCSI/unit)
  enable_testing()
endif()

#------------------------------------------------------------------------------#
# Add library.
#------------------------------------------------------------------------------#

add_library(Spec INTERFACE)
add_library(Spec::Spec ALIAS Spec)
target_include_directories(Spec
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
add_subdirectory(spec)

#------------------------------------------------------------------------------#
# Add apps.
#------------------------------------------------------------------------------#

add_subdirectory(app)
