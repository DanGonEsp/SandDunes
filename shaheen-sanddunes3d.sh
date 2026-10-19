#!/bin/bash -l

#SBATCH --job-name=Dune3d
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=8
#SBATCH --ntasks=1504
#SBATCH --ntasks-per-node=188
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/gonzald/SandDunes/Dune3d-%j.out
#SBATCH --error=/scratch/gonzald/SandDunes/Dune3d-%j.err
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1


UGSHELL=/project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell
WORKDIR=/scratch/gonzald/SandDunes
APPDIR=/project/k10105/gonzald/SandDunes/Files

cd "${APPDIR}" || exit 1

ARGS=(
    -ex Multiphase.lua
    -dir_name "${WORKDIR}"
    -dim 3
    -timeMethod limex
    -numRefs 5
    -numPreRefs 3
    -numTimeSteps 100
    -DT 3600.0
    -baseSolver SuperLu
    -riemman 3
)



# 1504 MPI ranks = 1 cases x 1504 spatial ranks
srun --exclusive --exact --ntasks=1504 --ntasks-per-node=188 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneRusanov3d-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneRusanov3d-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" &
PID1=$!


# 64 MPI ranks = 1 cases x 64 spatial ranks
#srun --exclusive --exact --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneEngquistOsher3d-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneEngquistOsher3d-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -riemman 3 &
#PID2=$!



STATUS=0
wait "${PID1}" || { echo "ERROR: 1504 MPI rank DuneRusanov3d experiment failed."; STATUS=1; }
#wait "${PID2}" ||  { echo "ERROR: 64  MPI rank DuneEngquistOsher experiment failed."; STATUS=1; }


if [ "${STATUS}" -ne 0 ]; then
    echo "=============================================="
    echo "One or more Dune 3d experiments failed."
    echo "Job ID: ${SLURM_JOB_ID}"
    echo "=============================================="
    exit 1
fi

echo "=============================================="
echo "All Dune 3d simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
