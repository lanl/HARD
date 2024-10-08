#SBATCH --job-name=multi

#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --partition=scaling

#SBATCH --output=res/output_%A_%a.out
#SBATCH --error=res/error_%A_%a.err
#SBATCH --time=03:00:00

export OMPI_MCA_btl=^openib
start_time=$(date +%s)
 
/usr/bin/time -v srun --cpu-bind=cores ../build/app/hard xcf_$SLURM_JOB_NUM_NODES.yaml

end_time=$(date +%s)

rm cf_$start_time.yaml

elapsed_time=$(($end_time - $start_time))
echo "Elapsed time: $elapsed_time seconds"
