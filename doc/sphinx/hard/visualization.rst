.. |br| raw:: html

   <br />

.. _flastro_visualization:

Visualization
*************

FlAstro has an integrated Catalyst
(https://www.paraview.org/Wiki/ParaView/Catalyst/Overview; https://docs.paraview.org/en/latest/Catalyst/index.html)
adaptor that creates output that can directly be opened in ParaView
(https://www.paraview.org/; https://docs.paraview.org/en/latest/UsersGuide/introduction.html).
**Note**: The Catalyst adaptor is only activated for 3D simulations. For 2D or 1D simulations,
FlAstro currently outputs raw text data files that can be read into gnuplot, for example. 

**How to build and run flastro with catalyst**

To build flastro with catalyst, set ``ENABLE_CATALYST=ON``.
Run flastro as you would w/o creating visualzation output, for example:

``mpirun -n 4 ./flastro -d 3 ../../configs/XYZ.yaml``

In the ``XYZ.yaml`` file you secify the physical dimension and resolution of the simulation box, the simulation time
or maximum number of steps, etc. To create output for paraview, you need to provide the path for the catalyst implementation
and specify a catalyst/paraview script, like ``gridwriter.py`` that is provided in the ``../../tools/`` directory along with
other example scripts that directly render pngs for different observables. To add further observables/fields to the pipeline
for visualization, see below.

3D visualizations with catalyst/paraview
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Example**: *Evolution of speed (magnitude of velocity) in a 3D Sod shock tube during the first 0.5 seconds.*

FlAstro's catalyst adaptor will put all visualization files into a directory ``datasets``.
It can be directly opened from your local paraview client. This is a visualization of an example run with
48 x 12 x12 cells, using ``gridwriter.py``:

.. image:: images/3d-flastro-sod-48x12x12cells.gif
  :width: 600
  :alt: Evolution of speed in a 3D Sod shock tube during the first 0.5 seconds.

2D visualizations using raw output and gnuplot/pm3d
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Flastro creates raw output files in ``tsv`` format that can be read into spead sheets or by 
tools like ``gnuplot``. Animations can be produced by plotting the data from each time step.

**Example**: Evolution of density in a 2D Sod shock tube during the first approx. 0.5 seconds.

.. image:: images/2d-flastro-sod-128cells.gif
  :width: 400
  :alt: Evolution of density in a 2D Sod shock tube during the first ~0.5 seconds.

**HowTo**:

#. Process all raw output files from a 2D-FlAstro run with ``gnuplot``:

   * the script ``make_png_from_raw_2D.sh`` in the ``../../tools/`` directory will create a png file from 
     every data file in a specified range

   * Example was made with: ``make_png_from_raw_2D.sh 1 200 128 4`` (for help, run script w/o arguments)

#. Combine pictures into an animation, for example using QuickTime or ImageMagick

   * QuickTime Player: load all image files ("Open image Sequence ...") and save as movie
   
   * ImageMagick (available on Darwin): ``convert *png 2d-flastro-sod-128cells.gif``


1D visualizations using raw output and gnuplot/pm3d
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create a plot of an observable at a certain time, you can directly plot columns from the raw output files.

**Example**: Density in a 1D Sod shock tube at time :math:`t = 0.2`
(https://en.wikipedia.org/wiki/Sod_shock_tube).

.. image:: images/output_1d_128_t0.2.png
  :width: 400
  :alt: Alternative text

To create a 3D plot of the time evolution of a scalar observable:

**Example**: Evolution of density in a 1D Sod shock tube during the first 0.2 seconds.

.. image:: images/output_1d_128_t0_0.2.png
  :width: 400
  :alt: Evolution of density in a 1D Sod shock tube during the first 0.2 seconds.

**HowTo**:

#. Concatenate all raw output files from a 1D-FlAstro run: ``cat output*raw > output_time_evolution.raw``
#. Run ``gnuplot``:

.. code-block::

  user@host:.../flastro/build/app$ gnuplot

  gnuplot> splot 'output_time_evolution.raw' u 1:($2/128):3 w pm3d at bs title "flastro 1D, 128 cells, sod evolution 0.0 < t < 0.2"
  gnuplot> set xlabel "time"
  gnuplot> set ylabel "x"
  gnuplot> set zlabel "density"

  [alternatively to all of the above:] 
  gnuplot> load "output_time_evolution.gnuplot"

  [use mouse to rotate to desired perspective]
  [to output png:]
  gnuplot> set term png
  gnuplot> set output "output_time_evolution.png"
  gnuplot> replot
  gnuplot> quit


How to add other fields/observables to the catalyst output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding new data to the catalyst/paraview pipeline requires editing all files in the catalyst module
(``catalyst/adaptor.[cc|hh]`` and ``catalyst/types[cc|hh]``), the catalyst tasks file (``tasks/catalyst.hh``),
and the files where you call catalyst tasks (``init.hh`` and ``analyze.hh``, in this example).

**Example: Adding the total energy**

**Note:** All necessary additions/changes for this example are available in the code repository, included in the commit
“How to add total energy to catalyst pipeline; documentation example”, and marked in the code with the ``<<< add total energy``
tag.

1. Physics fields/observables are referred to as “attributes” in the catalyst language.
So we have to create a new attribute and update it along with all the others. This is
done in ``catalyst/types.[hh|cc]`` and the ``update_attributes`` method in ``tasks/catalyst.hh``.
Note that total energy is a scalar field; in particular, it does not depend on the dimension.
(For vector fields, use ``velocity`` as an example accordingly).

``catalyst/types.hh``:

.. code-block::

  class catalyst_attributes 
  public:
    …
    void update_fields(double * cell_vector_vals1,
                      double * cell_scalar_vals1,
                      double * cell_scalar_vals2,
                      double * cell_scalar_vals3); // <<< add total energy
    …

  private:
    …
    std::vector<double> tot_energy; // <<< add total energy
    …
  };

``catalyst/types.cc``:

.. code-block::

  void catalyst_attributes::initialize(size_t number_of_points, size_t number_of_cells)
  {
    …
    this->tot_energy.resize(number_of_cells); // <<< add total energy
  }

  void catalyst_attributes::update_fields(double * cell_vector_vals1,
                                          double * cell_scalar_vals1,
                                          double * cell_scalar_vals2,
                                          double * cell_scalar_vals3) { // <<< add total energy
    for(size_t cell_i = 0; cell_i < num_cells_; cell_i++) {
      … 
      this->tot_energy[cell_i] = cell_scalar_vals3[cell_i]; // <<< add total energy
    }
    …
  } // end update_fields


``task/catalyst.hh``

.. code-block::

  template<std::size_t D>
  inline void update_attributes(single<catalyst_attributes>::accessor<rw> c_a,
    …
    field<double>::accessor<ro, ro> rE_a) // <<< add total energy
  {
    if constexpr(D == 3) {
      …
      auto rE = m.template mdcolex<is::cells>(rE_a); // <<< add total energy
      …
      std::vector<double> tot_energy_vals; // <<< add total energy

      for(auto k : m.template cells<ax::z, dm::quantities>()) {
       for(auto j : m.template cells<ax::y, dm::quantities>()) {
         for(auto i : m.template cells<ax::x, dm::quantities>()) {
            // update all scalar fields
            …
            tot_energy_vals.push_back(rE(i,j,k)); // <<< add total energy
            …
          } // for
        } // for
      } // for
      …
      c_a->update_fields(velocity_vals.data(),
                         density_vals.data(),
                         pressure_vals.data(),
                         tot_energy_vals.data()); // <<< add total energy
    }
    …
  }


Finally, pass the new field along when the ``update_attributes`` task is executed.

``init.hh`` and ``analyze.hh``

.. code-block::

  // Send initial problem state to catalyst
  execute<tasks::external::update_attributes<D>, mpi>(
    catalyst_data(pt), 0.0, s.m, s.r(s.m), s.u(s.m), s.p(s.m),
    s.rE(s.m)); // <<< add total energy

.. code-block::

  execute<tasks::external::update_attributes<D>, mpi>(
    catalyst_data(pt), cp.time(), s.m, s.r(s.m), s.u(s.m), s.p(s.m),
    s.rE(s.m)); // <<< add total energy

2. Now the new attribute has to be included in the conduit node that is created for catalyst/paraview.
This node is defined in ``catalyst/adaptor.cc``. Note that the total energy is “cell data”
(cells are called “elements” in catalyst language). We have to also add the new attribute to
the ``catalyst_attributes::execute`` method in ``catalyst/types.cc``.

``catalyst/adaptor.hh``

.. code-block::

  void execute(int cycle, double time, flastro::simple_cubic& grid,
    … 
    double* tot_energy_array); // <<< add total energy


``catalyst/adaptor.cc``

.. code-block::

  void catalyst_adaptor::execute(int cycle, double time, simple_cubic& grid,
    …
    double* tot_energy_array) // <<< add total energy
  {
  …
  // add total energy >>>
    // total energy is cell-data
    fields["total_energy/association"].set("element");
    fields["total_energy/topology"].set("mesh");
    fields["total_energy/volume_dependent"].set("false");
    fields["total_energy/values"].set_external(tot_energy_array, grid.get_number_of_cells());
  // <<< add total energy
  …
  }

``catalyst/types.cc``

.. code-block::

  void catalyst_attributes::execute(size_t step, double time, simple_cubic& lattice)
  {
    …
    catalyst_adaptor::execute(step, time, lattice, &this->velocity[0],
                                                  … 
                                                  &this->tot_energy[0]); // <<< add total energy
  }

Example for config file: 

.. code-block:: 
  
  catalyst:
    script: gridwriter.py
    implementation: paraview
    implementation_directory: /vast/home/thomasvogel/spack/var/spack/environments/vtk_viz/.spack-env/view/lib64/catalyst/



Feel free to contact Thomas (thomasvogel@lanl.gov) directly if you have any further questions.
