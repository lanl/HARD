#!/bin/bash

# Path to the executable
mpi_executable="$1"
hard_executable="$2"
config_file="$3"
additional_args="$4"

# Run the executable twice
rm *.csv
ls
echo "Starting.."
OMP_NUM_THREADS=1 "$mpi_executable" -np 1 "$hard_executable" "$config_file" $additional_args
echo "Mid..."
# rename before getting overwritten by next command
mv output-*-0-00002.csv r1.csv
ls
OMP_NUM_THREADS=1 "$mpi_executable" -np 16 "$hard_executable" "$config_file" $additional_args
ls

# Stitch together 16 csv outputs to 1 (mp.csv)
> mp.csv

for i in {0..15}; do
    # Exclude header and first 3 columns
    sed '1d;/^$/d' "output-1D-${i}-00002.csv"  | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $5, $6, $7}' >> mp.csv
done

# Create a temporary version of r1.csv excluding the header lines
sed '1d;/^$/d' r1.csv  | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $5, $6, $7}' > r1_temp.csv

diff mp.csv r1_temp.csv > /dev/null
exit_status=$?

if [ $exit_status -eq 0 ]; then
    echo "Files are the same."
else
    echo "Files are different."
    diff mp.csv r1_temp.csv
fi

# Clean up temporary files
rm *.csv

# Exit with the status from diff (0 for no differences, 1 for differences)
exit $exit_status

