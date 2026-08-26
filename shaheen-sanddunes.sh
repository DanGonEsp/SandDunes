#!/bin/bash
#SBATCH -N 2
#SBATCH --partition=workq
#SBATCH -J SandDunes
#SBATCH --output=SandDunes-output.txt
#SBATCH --error=SandDunes-error.txt
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10105
#SBATCH -t 23:59:00

#OpenMP settings:
export OMP_NUM_THREADS=1


#run the application:
srun --hint=nomultithread --ntasks=384 --ntasks-per-node=192 --ntasks-per-socket=96 --cpus-per-task=1 --ntasks-per-core=1 --mem-bind=v,local --cpu-bind=threads /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 3  -timeMethod limex -numTimeSteps 100
