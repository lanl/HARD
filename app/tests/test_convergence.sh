#!/bin/bash

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
start="$3"
stop="$4"
config_file="$5"
additional_args="$6"

rm -rf output*
for INDX in $(seq $start $stop)
do
    mkdir output_$INDX
    cd output_$INDX
    $mpi_executable -np 1 $hard_executable $config_file $additional_args -r $INDX
    cd ..
done
make convergence CONFIG=$config_file

rm -rf output*

# Exit with the status from diff (0 for no differences, 1 for differences)
exit $exit_status

