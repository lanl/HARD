#!/bin/bash
set -euo pipefail

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
start="$3"
stop="$4"
config_file="$5"
additional_args="$6"

dimension="${additional_args[0]: -1}"

rm -rf output*
for INDX in $(seq $start $stop)
do
    mkdir output_$INDX
    cd output_$INDX
    $mpi_executable -np 1 $hard_executable $config_file $additional_args -r $INDX
    cd ..
done
make convergence CONFIG=$config_file DIMENSION=$dimension

rm -rf output*

