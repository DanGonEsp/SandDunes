#!/bin/bash -l

#SBATCH --job-name=Dune2d
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/gonzald/SandDunes/Dune2d-%j.out
#SBATCH --error=/scratch/gonzald/SandDunes/Dune2d-%j.err
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
    -numRefs 7
    -numPreRefs 3
    -numTimeSteps 100
    -DT 3600.0

)



# 24 MPI ranks = 1 cases x 24 spatial ranks
srun --exclusive --exact --ntasks=24 --ntasks-per-node=24 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneRusanov2d-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneRusanov2d-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -riemman 1 &
PID1=$!


# 24 MPI ranks = 1 cases x 24 spatial ranks
srun --exclusive --exact --ntasks=24 --ntasks-per-node=24 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneEngquistOsher2d-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneEngquistOsher2d-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -riemman 3 &
PID2=$!



STATUS=0
wait "${PID1}" || { echo "ERROR: 24 MPI rank DuneRusanov2d experiment failed."; STATUS=1; }
wait "${PID2}" ||  { echo "ERROR: 24  MPI rank DuneEngquistOsher2d experiment failed."; STATUS=1; }


if [ "${STATUS}" -ne 0 ]; then
    echo "=============================================="
    echo "One or more Dune 2d experiments failed."
    echo "Job ID: ${SLURM_JOB_ID}"
    echo "=============================================="
    exit 1
fi

echo "=============================================="
echo "All Dune 2d simulations finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
