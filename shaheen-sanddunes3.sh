#!/bin/bash -l

#SBATCH --job-name=Avalanche3
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/gonzald/SandDunes/Avalanche3-%j.out
#SBATCH --error=/scratch/gonzald/SandDunes/Avalanche3-%j.err
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=ALL


export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

module load paraview/6.1.0-mesa

UGSHELL=/project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell
WORKDIR=/scratch/gonzald/SandDunes
APPDIR=/project/k10105/gonzald/SandDunes/Files

cd "${APPDIR}" || exit 1

ARGS=(
    -ex Avalanche.lua
    -dir_name "${WORKDIR}"
    -boolData false  #Always false#
    -dim 2
    -numProc 1
    -simCase 1
    -simCaseBnd 1
    -timeMethod limex
    -numTimeSteps 100
    -DT 1000.0
)

# ============================================================
# Weak-scaling simulations
# One simulation at a time
# ============================================================

srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numRefs 6 -numPreRefs 5

srun --exclusive --ntasks=16 --ntasks-per-node=16 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numRefs 5 -numPreRefs 4

srun --exclusive --ntasks=4 --ntasks-per-node=4 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numRefs 4 -numPreRefs 3

srun --exclusive --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numRefs 3 -numPreRefs 2


echo "=============================================="
echo "All Avalanche simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
