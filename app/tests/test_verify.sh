#!/bin/bash
set -euo pipefail

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
config_file="$3"
additional_args="$4"

$mpi_executable -np 1 $hard_executable $config_file $additional_args && \
  make verify CONFIG=$config_file
rm -rf *.raw

