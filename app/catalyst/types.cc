// Derived from FEDataStructures.h
//   SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
//   SPDX-License-Identifier: BSD-3-Clause
//   https://gitlab.kitware.com/paraview/paraview.git
//   https://gitlab.kitware.com/paraview/paraview/-/tree/master/Examples/Catalyst2/CxxFullExample

#include "types.hh"
#include "../types.hh"
#include "adaptor.hh"

#include <iostream>
#include <iterator>
#include <mpi.h>

#include <flecsi/flog.hh>

using namespace flecsi;
using namespace hard;

simple_cubic::simple_cubic() = default;

void
simple_cubic::initialize(const unsigned int num_points[3],
  const double spacing[3]) {
  if(num_points[0] == 0 || num_points[1] == 0 || num_points[2] == 0) {
    std::cerr
      << "Must have a non-zero amount of lattice points in each direction.\n";
  }

  unsigned int start_x_point =
    color_cell_id_begin[0]; // vertex id at beginning of cell = cell id
  unsigned int end_x_point =
    color_cell_id_end[0] + 1; // vertex id at end of cell = cell_id + 1
  unsigned int start_y_point = color_cell_id_begin[1];
  unsigned int end_y_point = color_cell_id_end[1] + 1;
  unsigned int start_z_point = color_cell_id_begin[2];
  unsigned int end_z_point = color_cell_id_end[2] + 1;

  /*   std::cout << "simple_cubic::initialize(), rank " << my_rank << ":"
               << " x: " << start_x_point << "-" << end_x_point
               << " y: " << start_y_point << "-" << end_y_point
               << " z: " << start_z_point << "-" << end_z_point << std::endl;
   */

  // create the points -- slowest in the z and fastest in the x directions
  double coord[3] = {0, 0, 0};
  int point_count = 0;
  for(unsigned int k = start_z_point; k <= end_z_point; k++) {
    coord[2] = k * spacing[2];
    for(unsigned int j = start_y_point; j <= end_y_point; j++) {
      coord[1] = j * spacing[1];
      for(unsigned int i = start_x_point; i <= end_x_point; i++) {
        coord[0] = i * spacing[0];
        // add the coordinate to the end of the vector
        std::copy(coord, coord + 3, std::back_inserter(this->points));
        point_count++;
      }
    }
  }

  // create the hex cells
  unsigned int cell_points[8];
  unsigned int num_x_points = end_x_point - start_x_point + 1;
  unsigned int num_y_points = end_y_point - start_y_point + 1;
  unsigned int num_z_points = end_z_point - start_z_point + 1;

  flog(info) << "simple_cubic::initialize() - create hex cells" << std::endl;
  // every rank has its own cell array starting at cell[0] with point ids 0, 1,
  // etc; just the physical _coordinates_ of those points are different in each
  // rank
  for(unsigned int k = 0; k < num_z_points - 1; k++) {
    for(unsigned int j = 0; j < num_y_points - 1; j++) {
      for(unsigned int i = 0; i < num_x_points - 1; i++) {
        /*
        See
        https://llnl-conduit.readthedocs.io/en/v0.7.1/blueprint_mesh.html#hexs
        for vertex connectivity: (0,0,0) -> (1,0,0) -> (1,1,0) -> (0,1,0) ->
                                 (0,0,1) -> (1,0,1) -> (1,1,1) -> (0,1,1)
        */
        cell_points[0] = k * num_y_points * num_x_points + j * num_x_points + i;
        cell_points[1] =
          k * num_y_points * num_x_points + j * num_x_points + i + 1;
        cell_points[2] =
          k * num_y_points * num_x_points + (j + 1) * num_x_points + i + 1;
        cell_points[3] =
          k * num_y_points * num_x_points + (j + 1) * num_x_points + i;
        cell_points[4] =
          (k + 1) * num_y_points * num_x_points + j * num_x_points + i;
        cell_points[5] =
          (k + 1) * num_y_points * num_x_points + j * num_x_points + i + 1;
        cell_points[6] = (k + 1) * num_y_points * num_x_points +
                         (j + 1) * num_x_points + i + 1;
        cell_points[7] =
          (k + 1) * num_y_points * num_x_points + (j + 1) * num_x_points + i;

        std::copy(
          cell_points, cell_points + 8, std::back_inserter(this->cells));
      }
    }
  }
}

size_t
simple_cubic::get_number_of_points() {
  return this->points.size() / 3;
}

size_t
simple_cubic::get_number_of_cells() {
  return this->cells.size() / 8;
}

double *
simple_cubic::get_points_array() {
  if(this->points.empty()) {
    return nullptr;
  }
  return &(this->points[0]);
}

double *
simple_cubic::get_point(size_t point_id) {
  if(point_id >= this->points.size()) {
    return nullptr;
  }
  return &(this->points[point_id * 3]);
}

unsigned int *
simple_cubic::get_cell_points(size_t cell_id) {
  if(cell_id >= this->cells.size()) {
    return nullptr;
  }
  return &(this->cells[cell_id * 8]);
}

catalyst_attributes::catalyst_attributes() {
  this->grid_ptr = nullptr;
}

void
catalyst_attributes::initialize(size_t number_of_points,
  size_t number_of_cells) {
  num_cells_ = number_of_cells;

  flog(info) << "Initialize velocity for catalyst with " << number_of_points
             << " elements" << std::endl;
  flecsi::flog::flush();

  this->velocity.resize(number_of_cells * 3);
  this->density.resize(number_of_cells);
  this->pressure.resize(number_of_cells);
  this->tot_energy.resize(number_of_cells);
  this->rad_energy.resize(number_of_cells);
  // <<< add variables here for catalyst
}

void
catalyst_attributes::update_fields(double * cell_vector_vals1,
  double * cell_scalar_vals1,
  double * cell_scalar_vals2,
  double * cell_scalar_vals3,
  double * cell_scalar_vals4) { // <<< add variables here for catalyst
  for(size_t cell_i = 0; cell_i < num_cells_; cell_i++) {
    this->density[cell_i] = cell_scalar_vals1[cell_i];
    this->pressure[cell_i] = cell_scalar_vals2[cell_i];
    this->tot_energy[cell_i] = cell_scalar_vals3[cell_i];
    this->rad_energy[cell_i] = cell_scalar_vals4[cell_i];
    // <<< add variables here for catalyst
  }
  for(size_t i = 0; i < (3 * num_cells_); i++) {
    this->velocity[i] = cell_vector_vals1[i];
  }
} // end update_fields

void
catalyst_attributes::execute(size_t step, double time, simple_cubic & lattice) {
  flog(info) << "catalyst_attributes::execute(): Send data to catalyst."
             << std::endl;
  catalyst_adaptor::execute(step,
    time,
    lattice,
    &this->velocity[0],
    &this->density[0],
    &this->pressure[0],
    &this->tot_energy[0],
    &this->rad_energy[0]); // <<< add variables here for catalyst
} // end execute

void
catalyst_attributes::finalize() {
  // here one could do things like releasing allocated memory,
  // performing certain analyses of the final state of the physical fields, etc.
}

double *
catalyst_attributes::get_velocity_array() {
  if(this->velocity.empty()) {
    return nullptr;
  }
  return &this->velocity[0];
}

double *
catalyst_attributes::get_pressure_array() {
  if(this->pressure.empty()) {
    return nullptr;
  }
  return &this->pressure[0];
}
