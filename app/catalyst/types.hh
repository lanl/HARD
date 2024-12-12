// Derived from FEDataStructures.h
//   SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
//   SPDX-License-Identifier: BSD-3-Clause
//   https://gitlab.kitware.com/paraview/paraview.git
//   https://gitlab.kitware.com/paraview/paraview/-/tree/master/Examples/Catalyst2/CxxFullExample

#ifndef FLASTRO_CATALYST_TYPES_HH
#define FLASTRO_CATALYST_TYPES_HH

#include <cstddef>
#include <vector>

namespace flastro {

class simple_cubic
{
public:
  simple_cubic();
  void initialize(const unsigned int num_points[3], const double spacing[3]);
  size_t get_number_of_points();
  size_t get_number_of_cells();
  double * get_points_array();
  double * get_point(size_t point_id);
  unsigned int * get_cell_points(size_t cell_id);

private:
  std::vector<double> points; // physical coordinates of lattice points
  std::vector<unsigned int>
    cells; // stores lattice point/vertex IDs that define each cell
};

class catalyst_attributes
// call it "attributes" to have consistent language with catalyst
// a "mesh" is a "grid" and physical fields (pressure, velocity, density, ...)
// are "attributes"
{
  // A class for generating and storing point and cell fields.
  // Velocity is (at the moment) stored at the points and pressure is stored
  // for the cells.
public:
  catalyst_attributes();
  void initialize(size_t number_of_points, size_t number_of_cells);
  void update_fields(double * cell_vector_vals1,
    double * cell_scalar_vals1,
    double * cell_scalar_vals2,
    double * cell_scalar_vals3,
    double * cell_scalar_vals4); // <<< add variables here for catalyst
  void execute(size_t step, double time, simple_cubic & lattice);
  void finalize();
  double * get_velocity_array();
  double * get_pressure_array();

private:
  std::size_t num_cells_;
  std::vector<double> velocity;
  std::vector<double> density;
  std::vector<double> pressure;
  std::vector<double> tot_energy;
  std::vector<double> rad_energy; // <<< add variables here for catalyst

  simple_cubic * grid_ptr;
};

inline simple_cubic lattice;
inline unsigned int number_vertices[3];
inline double lattice_spacing[3] = {1.0, 1.0, 1.0};
inline unsigned int
  number_colors[3]; // stores number of colors/blocks in each direction (x,y,z)
inline unsigned int
  color_coordinate[3]; // x,y,z coordinates of color blocks as a whole
inline unsigned int color_cell_id_begin[3]; // global id of cell at the
                                            // beginning of color block in x,y,z
inline unsigned int
  color_cell_id_end[3]; // global id of cell at the end of color block in x,y,z
inline unsigned int
  my_rank; // to mimic MPI language, stores flecsi process id / "rank"

} // namespace flastro
#endif // FLASTRO_CATALYST_TYPES_HH
