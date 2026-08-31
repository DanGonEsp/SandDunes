#!/bin/bash -l

#SBATCH --job-name=AvalancheSmall
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=1
#SBATCH --ntasks=124
#SBATCH --ntasks-per-node=124
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

UGSHELL=/project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell
APPDIR=/project/k10105/gonzald/SandDunes/Files
WORKDIR=/scratch/gonzald/SandDunes

cd "${APPDIR}" || exit 1

COMMON_ARGS=(
    -ex Avalanche.lua
    -dim 2
    -simCaseBnd 1
    -timeMethod limex
    -numProc 4
    -numRefs 4
    -numPreRefs 3
    -numTimeSteps 100
    -DT 10.0
)

# 64 MPI ranks = 4 cases x 16 spatial ranks
srun --exclusive --exact --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --output="${WORKDIR}/Avalanche64-${SLURM_JOB_ID}.out" --error="${WORKDIR}/Avalanche64-${SLURM_JOB_ID}.err" "${UGSHELL}" "${COMMON_ARGS[@]}" -dir_name "${WORKDIR}/Scaling_64" &
PID64=$!

# 32 MPI ranks = 4 cases x 8 spatial ranks
srun --exclusive --exact --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --output="${WORKDIR}/Avalanche32-${SLURM_JOB_ID}.out" --error="${WORKDIR}/Avalanche32-${SLURM_JOB_ID}.err" "${UGSHELL}" "${COMMON_ARGS[@]}" -dir_name "${WORKDIR}/Scaling_32" &
PID32=$!

# 16 MPI ranks = 4 cases x 4 spatial ranks
srun --exclusive --exact --ntasks=16 --ntasks-per-node=16 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --output="${WORKDIR}/Avalanche16-${SLURM_JOB_ID}.out" --error="${WORKDIR}/Avalanche16-${SLURM_JOB_ID}.err" "${UGSHELL}" "${COMMON_ARGS[@]}" -dir_name "${WORKDIR}/Scaling_16" &
PID16=$!

# 8 MPI ranks = 4 cases x 2 spatial ranks
srun --exclusive --exact --ntasks=8 --ntasks-per-node=8 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --output="${WORKDIR}/Avalanche8-${SLURM_JOB_ID}.out" --error="${WORKDIR}/Avalanche8-${SLURM_JOB_ID}.err" "${UGSHELL}" "${COMMON_ARGS[@]}" -dir_name "${WORKDIR}/Scaling_8" &
PID8=$!

# 4 MPI ranks = 4 cases x 1 spatial rank
srun --exclusive --exact --ntasks=4 --ntasks-per-node=4 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --output="${WORKDIR}/Avalanche4-${SLURM_JOB_ID}.out" --error="${WORKDIR}/Avalanche4-${SLURM_JOB_ID}.err" "${UGSHELL}" "${COMMON_ARGS[@]}" -dir_name "${WORKDIR}/Scaling_4" &
PID4=$!

STATUS=0

wait "${PID64}" || { echo "ERROR: 64 MPI rank experiment failed."; STATUS=1; }
wait "${PID32}" || { echo "ERROR: 32 MPI rank experiment failed."; STATUS=1; }
wait "${PID16}" || { echo "ERROR: 16 MPI rank experiment failed."; STATUS=1; }
wait "${PID8}"  || { echo "ERROR: 8 MPI rank experiment failed."; STATUS=1; }
wait "${PID4}"  || { echo "ERROR: 4 MPI rank experiment failed."; STATUS=1; }

if [ "${STATUS}" -ne 0 ]; then
    echo "=============================================="
    echo "One or more AvalancheSmall experiments failed."
    echo "Job ID: ${SLURM_JOB_ID}"
    echo "=============================================="
    exit 1
fi

echo "=============================================="
echo "All AvalancheSmall scaling experiments finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
