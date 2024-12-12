// Derived from CatalystAdaptor.h
//   SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
//   SPDX-License-Identifier: BSD-3-Clause
//   https://gitlab.kitware.com/paraview/paraview.git
//   https://gitlab.kitware.com/paraview/paraview/-/tree/master/Examples/Catalyst2/CxxFullExample

#ifndef FLASTRO_CATALYST_ADAPTOR_HH
#define FLASTRO_CATALYST_ADAPTOR_HH

#include "types.hh"

#include <catalyst.hpp>

#include <iostream>
#include <string>

/*Thomas Vogel:
I am separating the original CatalystAdaptor.h from Kitware (that contains
declarations and definitions) into a .h file (forward declarations) and a .cxx
file (definitions) to avoid "multiple definitions" error during linking.
Definitions are needed in multiple source files.
*/

/**
 * The namespace hold wrappers for the three main functions of the catalyst API
 * - initialize
 * - execute
 * - finalize
 * Although not required it often helps with regards to complexity to collect
 * catalyst calls under a class /namespace.
 */

namespace catalyst_adaptor {
// forward declarations
void initialize();
void execute(int cycle,
  double time,
  flastro::simple_cubic & grid,
  double * velocity_array,
  double * density_array,
  double * pressure_array,
  double * tot_energy_array,
  double * rad_energy_array); // <<< add variables here for catalyst
void finalize();
} // namespace catalyst_adaptor

#endif // FLASTRO_CATALYST_ADAPTOR_HH
