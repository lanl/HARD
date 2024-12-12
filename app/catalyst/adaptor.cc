// Derived from CatalystAdaptor.cxx
//   SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
//   SPDX-License-Identifier: BSD-3-Clause
//   https://gitlab.kitware.com/paraview/paraview.git
//   https://gitlab.kitware.com/paraview/paraview/-/tree/master/Examples/Catalyst2/CxxFullExample

#include "adaptor.hh"
#include "../options.hh"

#include <yaml-cpp/yaml.h>

using namespace flastro;

/*Thomas Vogel:
I am separating the original CatalystAdaptor.h from Kitware (that contains
declarations and definitions) apart into a .h file (forward declarations) and a
.cxx file (definitions) to avoid "multiple definitions" error during linking.
Definitions are needed in multiple source files.
*/

/**
 * In this example, we show how we can use Catalysts's C++
 * wrapper around conduit's C API to create Conduit nodes.
 * This is not required. A C++ adaptor can just as
 * conveniently use the Conduit C API to setup the
 * `conduit_node`. However, this example shows that one can
 * indeed use Catalyst's C++ API, if the developer so chooses.
 */
void
catalyst_adaptor::initialize() {
  YAML::Node config = YAML::LoadFile(opt::config.value());

  // Populate the catalyst_initialize argument based on the "initialize"
  // protocol [1]. [1]
  // https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-initialize
  conduit_cpp::Node node;

  // Using the arguments given to the driver set the filename for the catalyst
  // script and pass the rest of the arguments as arguments of the script
  // itself. To retrieve these  arguments from the script  use the `get_args()`
  // method of the paraview catalyst module [2]
  // [2]
  // https://kitware.github.io/paraview-docs/latest/python/paraview.catalyst.html
  node["catalyst/scripts/script/filename"].set_string(
    config["catalyst"]["script"].as<std::string>());

  // if the catalyst script comes after all other application arguments,
  // this can be done to pass arguments to the catalyst script
  // argn would be the position of the script
  /*
  for (int cc = argn+1; cc < argc; ++cc)
  {
    conduit_cpp::Node list_entry =
  node["catalyst/scripts/script/args"].append(); list_entry.set(argv[cc]);
  }
  */

  // For this example we hardcode the implementation name to "paraview" and
  // define the "PARAVIEW_IMPL_DIR" during compilation time (see the
  // accompanying CMakeLists.txt). We could however defined them via
  // environmental variables  see [1].
  node["catalyst_load/implementation"] =
    config["catalyst"]["implementation"].as<std::string>();
  node["catalyst_load/search_paths/paraview"] =
    config["catalyst"]["implementation_directory"].as<std::string>();
  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if(err != catalyst_status_ok) {
    std::cerr << "Failed to initialize Catalyst: " << err << std::endl;
  }
}

void
catalyst_adaptor::execute(int cycle,
  double time,
  simple_cubic & grid,
  double * velocity_array,
  double * density_array,
  double * pressure_array,
  double * tot_energy_array,
  double * rad_energy_array) { // <<< add variables here for catalyst
  // Populate the catalyst_execute argument based on the "execute" protocol [3].
  // [3]
  // https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-execute

  conduit_cpp::Node exec_params;

  // State: Information about the current iteration. All parameters are
  // optional for catalyst but downstream filters may need them to execute
  // correctly.

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(cycle);
  state["time"].set(time);
  state["multiblock"].set(1);

  // Channels: Named data-sources that link the data of the simulation to the
  // analysis pipeline in other words we map the simulation datastructures to
  // the ones expected by ParaView.  In this example we use the Mesh Blueprint
  // to describe data see also bellow.

  // Add channels.
  // We only have 1 channel here. Let's name it 'grid'.
  auto channel = exec_params["catalyst/channels/grid"];

  // Since this example is using Conduit Mesh Blueprint to define the mesh,
  // we set the channel's type to "mesh".
  channel["type"].set("mesh");

  // now create the mesh.
  auto mesh = channel["data"];

  // populate the data node following the Mesh Blueprint [4]
  // [4] https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html

  // start with coordsets (of course, the sequence is not important, just make
  // it easier to think in this order).
  mesh["coordsets/coords/type"].set("explicit");
  // .set_external passes just the pointer  to the analysis pipeline allowing
  // thus for zero-copy data conversion see
  // https://llnl-conduit.readthedocs.io/en/latest/tutorial_cpp_ownership.html
  mesh["coordsets/coords/values/x"].set_external(grid.get_points_array(),
    grid.get_number_of_points(),
    /*offset=*/0,
    /*stride=*/3 * sizeof(double));
  mesh["coordsets/coords/values/y"].set_external(grid.get_points_array(),
    grid.get_number_of_points(),
    /*offset=*/sizeof(double),
    /*stride=*/3 * sizeof(double));
  mesh["coordsets/coords/values/z"].set_external(grid.get_points_array(),
    grid.get_number_of_points(),
    /*offset=*/2 * sizeof(double),
    /*stride=*/3 * sizeof(double));

  // Next, add topology
  mesh["topologies/mesh/type"].set("unstructured");
  mesh["topologies/mesh/coordset"].set("coords");
  mesh["topologies/mesh/elements/shape"].set("hex");
  mesh["topologies/mesh/elements/connectivity"].set_external(
    grid.get_cell_points(0), grid.get_number_of_cells() * 8);

  // Finally, add fields.

  // First component of the path is the name of the field .
  // The rest are described in
  // https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#fields
  // under the Material-Independent Fields section.
  auto fields = mesh["fields"];

  // if a field lives on the vertices, you can do this:
  // fields["velocity/association"].set("vertex");
  // if a field is defined a the cell level, do this ("element" = "cell")
  fields["velocity/association"].set("element");
  fields["velocity/topology"].set("mesh");
  fields["velocity/volume_dependent"].set("false");

  // velocity is stored in non-interlaced form (unlike points).
  fields["velocity/values/x"].set_external(
    velocity_array, grid.get_number_of_cells(), /*offset=*/0);
  fields["velocity/values/y"].set_external(velocity_array,
    grid.get_number_of_cells(),
    /*offset=*/grid.get_number_of_cells() * sizeof(double));
  fields["velocity/values/z"].set_external(velocity_array,
    grid.get_number_of_cells(),
    /*offset=*/grid.get_number_of_cells() * sizeof(double) * 2);

  // density is cell-data.
  fields["density/association"].set("element");
  fields["density/topology"].set("mesh");
  fields["density/volume_dependent"].set("false");
  fields["density/values"].set_external(
    density_array, grid.get_number_of_cells());

  // pressure is cell-data.
  fields["pressure/association"].set("element");
  fields["pressure/topology"].set("mesh");
  fields["pressure/volume_dependent"].set("false");
  fields["pressure/values"].set_external(
    pressure_array, grid.get_number_of_cells());

  // total energy is cell-data
  fields["total_energy/association"].set("element");
  fields["total_energy/topology"].set("mesh");
  fields["total_energy/volume_dependent"].set("false");
  fields["total_energy/values"].set_external(
    tot_energy_array, grid.get_number_of_cells());

  // radiation energy is cell-data
  fields["rad_energy/association"].set("element");
  fields["rad_energy/topology"].set("mesh");
  fields["rad_energy/volume_dependent"].set("false");
  fields["rad_energy/values"].set_external(
    rad_energy_array, grid.get_number_of_cells());

  // <<< add variables here for catalyst

  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if(err != catalyst_status_ok) {
    std::cerr << "Failed to execute Catalyst: " << err << std::endl;
  }
}

// Although no arguments are passed for catalyst_finalize  it is required in
// order to release any resources the ParaViewCatalyst implementation has
// allocated.
void
catalyst_adaptor::finalize() {
  conduit_cpp::Node node;
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if(err != catalyst_status_ok) {
    std::cerr << "Failed to finalize Catalyst: " << err << std::endl;
  }
}
