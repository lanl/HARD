.. |br| raw:: html

   <br />

.. _hard_visualization:

Visualization
*************

HARD supports two output formats for visualization and analysis.

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

Using VTK Output with ParaView
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HARD creates one master file per timestep: ``output-{D}D-{timestep}.pvti``

To visualize:

1. Open ParaView
2. File → Open → select the ``.pvti`` file
3. Click "Apply"

The master file automatically loads data from all MPI ranks.

Configuration Example
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

  config = {
    "problem": "sod",
    "output_frequency": 10,
    "output_method": "vti",  # or "csv"
    ...
  }
