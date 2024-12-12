#!/bin/bash
#SBATCH --job-name=KH
#SBATCH --exclusive
#SBATCH -p scaling
#SBATCH -N 4
#SBATCH --ntasks-per-node=36
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.out
#SBATCH -t 03:00:00

export OMPI_MCA_btl=^openib
start_time=$(date +%s)
 
/usr/bin/time -v srun --cpu-bind=cores ./app/flastro \
 ../configs/radiative-kh-instability.yaml ../tools/gridwriter.py -d 3

# /usr/bin/time -v srun --cpu-bind=cores ./app/flastro \
#  ../configs/kh.yaml ../tools/gridwriter.py -d 3

end_time=$(date +%s)
rm cf_$start_time.yaml
elapsed_time=$(($end_time - $start_time))
echo "Elapsed time: $elapsed_time seconds"
