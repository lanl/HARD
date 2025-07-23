#!/bin/bash
set -euo pipefail

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
config_file="$3"
additional_args="$4"

dimension="${additional_args[0]: -1}"

rm -rf *.csv
$mpi_executable -np 1 $hard_executable $config_file $additional_args && \
  make verify CONFIG=$config_file DIMENSION=$dimension
rm -rf *.csv

if [ -d ./artifacts ]; then
    # Save the plots as artifacts, suppressing erros for
    # those cases that do not produce plots
    mv *.png ./artifacts 2>/dev/null
fi
