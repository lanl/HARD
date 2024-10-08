#!/bin/bash

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
config_file="$3"
additional_args="$4"

# Run the executable twice
rm *.raw
ls
echo "Starting.."
OMP_NUM_THREADS=1 "$mpi_executable" -np 1 "$hard_executable" "$config_file" $additional_args
echo "Mid..."
# rename before getting overwritten by next command
mv output-00002-*-0.raw r1.raw
ls
OMP_NUM_THREADS=1 "$mpi_executable" -np 16 "$hard_executable" "$config_file" $additional_args
ls

# Stitch together 16 raw outputs to 1 (mp.raw)
> mp.raw

for i in {0..15}; do
    # Exclude header and first 3 columns
    sed '1,5d;/^$/d' "output-00002-1D-${i}.raw"  | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $5, $6, $7}' >> mp.raw
done

# Create a temporary version of r1.raw excluding the header lines
sed '1,5d;/^$/d' r1.raw  | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $5, $6, $7}' > r1_temp.raw

diff mp.raw r1_temp.raw > /dev/null
exit_status=$?

if [ $exit_status -eq 0 ]; then
    echo "Files are the same."
else
    echo "Files are different."
    diff mp.raw r1_temp.raw
fi

# Clean up temporary files
#rm *.raw

# Exit with the status from diff (0 for no differences, 1 for differences)
exit $exit_status

