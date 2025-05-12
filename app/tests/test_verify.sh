#!/bin/bash

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
config_file="$3"
additional_args="$4"

$mpi_executable -np 1 $hard_executable $config_file $additional_args
make verify CONFIG=$config_file
rm -rf *.raw

# Exit with the status from diff (0 for no differences, 1 for differences)
exit $exit_status

