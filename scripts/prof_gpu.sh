#!/bin/bash

#SBATCH --job-name=n1g2               # Job name
#SBATCH --partition=gpu               # Partition name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=2                    # Total number of MPI tasks
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --exclusive                   # Exclusive use of the nodes
#SBATCH --account=j19_cdss_g          # Account name
#SBATCH --output=hard-%j.out       # Standard output and error log

# activate env
cd ~/gpu-env/systems/chicoma/
source activate.sh craympich-mpi-clang/
cd ~/hard/build_ch/ 

date
hostname

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1

# Print out environment information for each rank
srun --output=hard_rank_%j_%t.out ./select_gpu nsys profile --output=chic_prof_%p --trace=nvtx,cuda,osrt ./app/hard ~/CONFIGS/rk.yaml -h 1


# Display job information
scontrol show job $SLURM_JOB_ID
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode

date    
