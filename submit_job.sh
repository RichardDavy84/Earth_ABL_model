#!/bin/bash
#SBATCH --account=nn9878k
#SBATCH --job-name=ABL
#SBATCH --time=00:20:00
#SBATCH --qos=devel
##SBATCH --partition=bigmem
#SBATCH --ntasks=1 --cpus-per-task=32
##SBATCH --mem-per-cpu=32G
#SBATCH --output=slurm_out.%j.out
#SBATCH --error=slurm_err.%j.err

cd /cluster/projects/nn9878k/hregan/ABL/Earth_ABL_model_parallel/test_cpupertask/test_32_Fram/

export OMP_NUM_THREADS=32

# vmstat 60&
srun --cpus-per-task=32 ABL.x 

