#!/bin/bash -l

#SBATCH --job-name=Dune
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/gonzald/SandDunes/Dune-%j.out
#SBATCH --error=/scratch/gonzald/SandDunes/Dune-%j.err
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
    -dim 2
    -timeMethod limex
    -numRefs 6
    -numPreRefs 3
    -numTimeSteps 100
    -DT 3600.0
)



# 64 MPI ranks = 1 cases x 64 spatial ranks
srun --exclusive --exact --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneRusanov-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneRusanov-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -riemman 1 &
PID1=$!


# 64 MPI ranks = 1 cases x 64 spatial ranks
srun --exclusive --exact --ntasks=64 --ntasks-per-node=64 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneEngquistOsher-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneEngquistOsher-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -riemman 3 &
PID2=$!



STATUS=0
wait "${PID1}" || { echo "ERROR: 64 MPI rank DuneRusanov experiment failed."; STATUS=1; }
wait "${PID2}" ||  { echo "ERROR: 64  MPI rank DuneEngquistOsher experiment failed."; STATUS=1; }


if [ "${STATUS}" -ne 0 ]; then
    echo "=============================================="
    echo "One or more Dune experiments failed."
    echo "Job ID: ${SLURM_JOB_ID}"
    echo "=============================================="
    exit 1
fi

echo "=============================================="
echo "All Dune simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
