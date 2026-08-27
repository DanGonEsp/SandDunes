#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=workq
#SBATCH -J SandDunes
#SBATCH --output=/scratch/gonzald/SandDunes/SandDunes-output.txt
#SBATCH --error=/scratch/gonzald/SandDunes/SandDunes-error.txt
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10105
#SBATCH -t 01:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1


#run the application:
srun -n 192 --ntasks-per-node=192 /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 3  -timeMethod limex -numTimeSteps 1 -

#End of File
