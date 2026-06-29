# HARD: Hydrodynamics And Radiative Diffusion

HARD is a radiation-hydrodynamics solver suite for the study of astrophysical phenomena.

HARD is based on the FleCSI framework and implemented on top of FleCSI-SP (FleCSI Specialization project).

# Copyright
© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so (Copyright request O4795).

# Spack build

The easiest way to build HARD is to use *spack*.

Clone the spack repo (version 1.0 or later required) and initialize:
```
$ git clone git@github.com:spack/spack.git $HOME/.spack
$ cd $HOME/.spack
$ source $HOME/.spack/share/spack/setup-env.sh
```

You can automate your spack setup by adding something like this to your
`$HOME/.bashrc` file:
```
# Setup spack environment
[ -f $HOME/.spack/share/spack/setup-env.sh ] && \
. $HOME/.spack/share/spack/setup-env.sh

export SPACK_EDITOR=vi
```
Notice that this also sets the *SPACK_EDITOR* environment variable. This
is useful for selecting the editor program spack will use for various
interface operations (used below).

Next, create a spack *environment* for HARD and activate it:
```
$ spack env create hard
$ spacktivate hard
```
Once you are in the *hard* environment, you can specify the repository:
```
$ spack repo add /PATH-TO-HARD-CLONE/spack_repo/hard
```
You can find the compiler we loaded earlier using:
```
$ spack compiler find
```

You can see the different options available for *hard* by using:
```
$ spack info hard
```

Available variants include:
- `+cuda` - Enable CUDA support for GPU acceleration
- `+radiation` - Enable radiation physics (default: on)
- `+tests` - Enable unit tests
- `+verification` - Enable physics verification tests
- `+format` - Enable code formatting target

Add *hard* to the environment:
```
$ spack add hard
```
The next step is to *concretize* the
new environment. Concretization solves an optimization problem that
takes all of the package dependencies in an environment and reconciles
version and subdependency compatibility:
```
$ spack concretize -f
```
The `-f` operation forces spack to re-concretize, which is not necessary
for a new, unconcretized environment. However, it doesn't hurt, and it
is useful in most invocations to insure that spack is using the most
up-to-date config information.

To finish up the spack portion of the build, we just need to *install*
the packages:
```
$ spack install --only dependencies
```
This will fetch all of the *hard* dependencies and build them in your
local spack tree (under *$HOME/.spack/opt* to be precise). This step
will take some time.

# Build

Now that we have setup our spack environment, the rest of the build is
really easy. First, change directory into your *hard* clone and create
a build directory:
```
$ cd PATH-TO-HARD-CLONE
$ mkdir build
```
Next, change directory into the build directory and run *cmake*:
```
$ cd build
$ cmake ..
```
This will configure your build with the default settings for HARD.

Available CMake options (use `-DOPTION=ON/OFF`):
- `ENABLE_HDF5` - Enable HDF5-based I/O for XDMF output (default: ON, requires parallel HDF5)
- `ENABLE_FORMAT` - Enable code formatting target (default: OFF, requires LLVM 20)
- `ENABLE_DOCUMENTATION` - Enable documentation generation (default: OFF)
- `ENABLE_UNIT_TESTS` - Enable unit tests (default: OFF)
- `ENABLE_VERIFICATION` - Enable physics verification tests (default: OFF, requires unit tests)
- `HARD_BENCHMARK_MODE` - Disable I/O and add time measurement for benchmarking (default: OFF)
- `HARD_ENABLE_LEGION_TRACING` - Enable Legion tracing (default: OFF)
- `HARD_WRITE_CONTROL_INFO` - Output control model graph at startup (default: ON, requires FleCSI Graphviz support)

When cmake has completed, simply run make:
```
$ make
```

# Output Formats

HARD supports multiple output formats for visualization in ParaView:

- **CSV** - Text-based output (always available)
- **VTK** - Parallel XML VTK format (.pvti/.vti files, always available)
- **XDMF+HDF5** - Binary HDF5 data with XDMF metadata (requires `ENABLE_HDF5=ON` and parallel HDF5)

For large-scale simulations, XDMF+HDF5 provides the most efficient I/O with full parallel write support.

# Running Simulations

HARD uses Python-based configuration files for problem setup. Example configuration files are provided in the `configs/` directory:

```
$ ./apps/hydro/hydro-{backend} -i configs/sod.py
```

Configuration files define initial conditions, boundary conditions, equation of state parameters, and output settings using Python syntax (replacing the previous YAML format).

# Code Formatting and Style

HARD enforces consistent code formatting and naming conventions to maintain code quality and readability.

## Formatting Requirements

All code must be formatted using `clang-format` before submitting a merge request. This is **mandatory** for any MR to be merged.

To format your code, run:
```
$ make format
```

This command will automatically format all source files according to the project's style guidelines.

## Naming Conventions

HARD uses **snake_case** for all variable, function, and file names. This convention must be followed consistently throughout the codebase.

Examples:
- Variables: `mass_density`, `velocity_field`, `time_step`
- Functions: `compute_flux()`, `update_state()`, `initialize_mesh()`
- Files: `riemann_solvers.hh`, `interface_fluxes.cc`

**Note:** Before submitting any merge request, ensure that:
1. Your code is properly formatted with `make format`
2. All naming follows the snake_case convention

<!-- vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 : -->
