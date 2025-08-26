.. |br| raw:: html

   <br />

.. _hard_input_files:

HARD YAML Input Files
************************************

HARD uses YAML files as input to describe the different problems and their configurations.
These files can be found in the `configs` directory at the root.
A subdirectory `ci_configs` is used to run the regression tests in the CI and should not be altered.

Unless specified, the units are in CGS.


Here are details on each parts of the yaml files:

`problem`
  This is the string that is picked up by the `init.hh` file to configure the problem. For now these options are available:
    * `sod`
    * `rankine-hugoniot`
    * `leblanc`
    * `acoustic-wave`
    * `kh-test`
    * `heating_and_cooling`
    * `sedov`
    * `implosion`
    * `rad-rh`
    * `lw-implosion`

Physics quantities
  Input that represent physical quantities: `gamma`, `kappa`, `mean_molecular_weight`

`eos`
  Define which EOS to be used from Singularity EOS
    * `ideal`
    * `spiner` (tabulated EOS, not tested)
    * `gruneisen` (not tested)

time and steps
  These variables handle the time and timestep of the simulation:
    * `t0`, start time for the simulation
    * `tf`, final time for the simulation (not supported for Legion backend)
    * `cfl`, Courant–Friedrichs–Lewy coefficient
    * `max_dt`, Maximum dt to be used in the simulation
    * `max_steps`, Maximum number of iterations to run the simulation for

logging and output
  * `log_frequency`, Set the frequency of iterations outputted on the terminal
  * `output_frequency`, Set the frequency of the output files

`color_distribution`
  An array representing how many colors per dimensions to assign.

`levels`
   This is an array representing the grid size as a power of two per dimension.
   As an example [8,8,4] would mean a grid of 2^8 * 2^8 * 2^4.
   If the code runs with the -d command line option, the extra dimensions are ignored.
`lowest_level`
  This defines the lowest level of resolution for the Geometric Multigrid.
  All dimensions decrease levels at the same rate, so the lowest resolution corresponds to whichever has one dimension arriving to the lowest level first
`coords`
  Coordinates of the lowest and highest points of the domain.
  We set two array per dimensions as, for example, a square:
  - [0.0, 0.0, 0.0]
  - [1.0, 1.0, 1.0]

`boundaries`
  The boundaries are defined per dimensions `xlow`, `xhigh`, `ylow`, `yhigh`, `zlow`, and `zhigh` For each several options are possible:
    * `periodic`
    * `flow`
    * `reflecting`
    * `dirichlet`

`catalyst`
  All these options are needed under catalyst to specify the location of the script and the paraview library.
    * `script`, path to the python script that interface the output data and catalyst. An example is provided in `tools/gridwriter.py`
    * `implementation`, only value for now is `paraview`
    * `implementation_directory`, path to the catalyst library compiled using the implementation. Something like `libcatalyst_paraview.so`


Some problems have specific input information. These are all defined under the `problem_parameters` group.

  * `acoustic_wave`:
    
    * `problem_parameters`, a list containing:
        * `r0`, the equilibrium density
        * `p0`, the equilibrium pressure
        * `amplitude`, the perturbation amplitudes
        * `scale`, the inverse of the wave number
  * `heating_and_cooling`:
    
    * `problem_parameters`, a list containing:
        * `fluid_mass_density`
        * `fluid_temperature`
        * `radiation_temperature`
  * `implosion`:
    
    * `problem_parameters`, a list containing:
        * `fluid_mass_density`
        * `fluid_temperature`
        * `radiation_temperature`
  * `sedov`:
    
    * `problem_parameters`, a list containing:
        * `hotspot_position`, an array of the position of the hotspot
        * `hotspot_radius`, the radius of the hotspot
        * `E_0`, the initial energy injected at the hotspot in Ergs

The other tests have setting hardcoded in their initialization files.
They can be found in `app/tasks/initial_data/`.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
