# Darwin

Get a `scaling` node on Darwin:

```
$ salloc -p scaling -t 01:00:00
```

Follow the instructions on the main README.md.
The following packages should be loaded in the environment before concretizing:

You can load the following compilers/mpi, and load them in Spack:

```
$ module load openmpi/5.0.2-gcc_13.2.0
$ spack compiler find
```

You can add the openmpi by adding it in the `package.py` file:

``` shell
$ vim ~/.spack/packages.py
```
And add:

``` yaml
packages:
  openmpi:
    externals:
    - spec: openmpi@5.0.2
      prefix: /projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0/
```

You can then run the concretizer.

## Format + Python

In order to use the python script and the clang format I would recommand loading the packages in the following order:

```
module load clang/13.0.0
spack env activate hard
module load miniconda3
```

# Chicoma
