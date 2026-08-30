#!/bin/bash -l

#SBATCH --job-name=Avalanche128
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


# ============================================================
# OpenMP settings
# ============================================================

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1


# ============================================================
# Application
# ============================================================

UGSHELL=/project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell
WORKDIR=/scratch/gonzald/SandDunes


# ============================================================
# Common Avalanche arguments
# ============================================================

ARGS=(
    -ex Avalanche.lua
    -dir_name "${WORKDIR}"
    -dim 2
    -simCaseBnd 1
    -timeMethod limex
    -numProc 1
    -numRefs 4
    -numPreRefs 3
    -numTimeSteps 100
    -DT 10.0
)



# ============================================================
# 128 MPI ranks
#
# numProc = 1
# No space-time splitting.
# One spatial communicator with 128 ranks.
# ============================================================


# 128 MPI ranks - simCase 1
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 1


# 128 MPI ranks - simCase 2
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 2


# 128 MPI ranks - simCase 3
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 3


# 128 MPI ranks - simCase 4
srun --exclusive --ntasks=128 --ntasks-per-node=128 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 4


# ============================================================
# 64 MPI ranks
#
# numProc = 1
# No space-time splitting.
# One spatial communicator with 64 ranks.
# ============================================================


# 64 MPI ranks - simCase 1
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 1


# 64 MPI ranks - simCase 2
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 2


# 64 MPI ranks - simCase 3
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 3


# 64 MPI ranks - simCase 4
srun --exclusive --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "${UGSHELL}" "${ARGS[@]}" -simCase 4


# ============================================================
# End
# ============================================================

echo "=============================================="
echo "All Avalanche simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="



#End of File
