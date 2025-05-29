# HARD: Hydrodynamics And Radiative Diffusion

HARD is a radiation-hydrodynamics solver suite for the study of astrophysical phenomena.

HARD is based on the FleCSI framework and implemented on top of FleCSI-SP (FleCSI Specialization project).

# Copyright
© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so (Copyright request O4795).

# Spack

The easiest way to build HARD is to use *spack*:

Clone the spack repo and initialize:
```
$ git clone git@github.com:spack/spack.git $HOME/.spack
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
$ spack repo add /PATH-TO-HARD-CLONE/spack-repo
```
You can see the different options available for *hard* by using:
```
$ spack info hard
```
And add *hard* to the environment:
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
When cmake has completed, simply run make:
```
$ make
```

<!-- vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 : -->
