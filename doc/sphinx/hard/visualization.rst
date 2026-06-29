.. |br| raw:: html

   <br />

.. _hard_visualization:

Visualization
*************

HARD supports three output formats for visualization and analysis.

Output Formats
~~~~~~~~~~~~~~

**CSV Output** (default)

Text-based output compatible with spreadsheet tools, gnuplot, and custom scripts.

.. code-block:: python

  "output_method": "csv",

**VTK Output**

Binary output for interactive 3D visualization in ParaView (https://www.paraview.org/).
Uses parallel VTK ImageData format (.vti/.pvti).

.. code-block:: python

  "output_method": "vti",

**XDMF Output** (requires HDF5)

Efficient binary output for large-scale simulations using XDMF+HDF5 format.
One HDF5 file per timestep with lightweight XDMF metadata.

.. code-block:: python

  "output_method": "xdmf",

Using VTK Output with ParaView
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HARD creates one master file per timestep: ``output-{D}D-{timestep}.pvti``

To visualize:

1. Open ParaView
2. File → Open → select the ``.pvti`` file
3. Click "Apply"

The master file automatically loads data from all MPI ranks.

Using XDMF Output with ParaView
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HARD creates two files per timestep:

- ``output-{D}D-{timestep}.h5`` - HDF5 binary data (all ranks write collectively)
- ``output-{D}D-{timestep}.xmf`` - XDMF metadata (lightweight XML)

To visualize:

1. Open ParaView
2. File → Open → select the ``.xmf`` file
3. Click "Apply"

ParaView reads metadata from the XDMF file and loads data from HDF5 on-demand.

**Benefits of XDMF:**

- Efficient binary storage with HDF5
- Parallel I/O (one file per timestep for all ranks)
- Smaller file sizes compared to VTI for large simulations
- Native ParaView support

Configuration Example
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

  config = {
    "problem": "sod",
    "output_frequency": 10,
    "output_method": "xdmf",  # or "vti" or "csv"
    ...
  }
