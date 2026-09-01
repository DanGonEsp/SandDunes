#!/bin/bash -l

#SBATCH --job-name=Avalanche2
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/gonzald/SandDunes/SandDunes-%j.out
#SBATCH --error=/scratch/gonzald/SandDunes/SandDunes-%j.err
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
    -simCaseBnd 1
    -timeMethod limex
    -numRefs 4
    -numPreRefs 3
    -numTimeSteps 100
    -DT 10.0
)

# ============================================================
# 128 MPI ranks
# One simulation at a time, using all 128 spatial ranks
# ============================================================

srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 1
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 2
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 3
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 4

# ============================================================
# 64 MPI ranks
# One simulation at a time, using 64 spatial ranks
# ============================================================

srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 1
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 2
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 3
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 1 -simCase 4


# ============================================================
# 128 MPI ranks
# One simulation at a time, using 32 spatial ranks
# ============================================================
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -numProc 4




echo "=============================================="
echo "All Avalanche simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
