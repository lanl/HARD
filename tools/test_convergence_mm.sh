#!/bin/bash
set -euo pipefail

# Path to the executable
mpi_executable="$1" 
hard_executable="$2"
start="$3"
stop="$4"
config_file="$5"
additional_args="$6"

dimension="$(printf '%s\n' "$additional_args" | sed -E 's/.*-d([123]).*/\1/')"

cleanup() {
    rm -rf output*
}
trap cleanup EXIT

rm -rf output*
rm -f contact_move_order_1_convergence-1D_*

for INDX in $(seq "$start" "$stop")
do
    mkdir "output_$INDX"
    cd "output_$INDX"
    "$mpi_executable" -n 1 "$hard_executable" "$config_file" "$additional_args" -r "$INDX"
    cd ..
done

set +e
python3 ../tools/convergence_verify_mm.py "$config_file" "$dimension" --plot
verify_status=$?
set -e

exit "$verify_status"
