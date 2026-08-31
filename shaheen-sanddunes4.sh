#!/bin/bash -l

#SBATCH --job-name=Avalanche4
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/gonzald/SandDunes/Avalanche4-%j.out
#SBATCH --error=/scratch/gonzald/SandDunes/Avalanche4-%j.err
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=ALL

set -e

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
    -file_name WeakSol
    -folder_name WeakSol
    -dim 2
    -numProc 1
    -simCase 1
    -simCaseBnd 1
    -timeMethod limex
    -numTimeSteps 100
    -DT 10.0
)

# ============================================================
# Weak-scaling simulation
# One simulation using 256 spatial MPI ranks
# 128 ranks per node across 2 nodes
# ============================================================

srun --exclusive --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numRefs 7 -numPreRefs 3

echo "=============================================="
echo "All Avalanche simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
