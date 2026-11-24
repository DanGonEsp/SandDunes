#!/bin/bash -l

#SBATCH --job-name=Dune2d
#SBATCH --partition=workq
#SBATCH --account=k10215
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=23:59:59
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
    -numTimeSteps 1000
    -outputFactor 10
    -DT 3600.0
    -DTLimex 1e-04
)




srun --exclusive --exact --ntasks=8 --ntasks-per-node=8 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneRusanov_01-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneRusanov_01-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -pre_smoothLim 1 -post_smoothLim 1 -riemman 1 &
PID1=$!


srun --exclusive --exact --ntasks=16 --ntasks-per-node=16 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneRusanov_02-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneRusanov_02-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -pre_smoothLim 2 -post_smoothLim 2 -riemman 1 &
PID2=$!

srun --exclusive --exact --ntasks=24 --ntasks-per-node=24 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneRusanov_03-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneRusanov_03-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -pre_smoothLim 2 -post_smoothLim 2 -riemman 1 &
PID3=$!


srun --exclusive --exact --ntasks=8 --ntasks-per-node=8 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneEngquistOsher2d_01-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneEngquistOsher2d_01-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -pre_smoothLim 1 -post_smoothLim 1 -riemman 3 &
PID4=$!


srun --exclusive --exact --ntasks=16 --ntasks-per-node=16 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneEngquistOsher2d_02-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneEngquistOsher2d_02-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -pre_smoothLim 2 -post_smoothLim 2 -riemman 3 &
PID5=$!

srun --exclusive --exact --ntasks=24 --ntasks-per-node=24 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores --mem=0 --output="${WORKDIR}/DuneEngquistOsher2d_03-${SLURM_JOB_ID}.out" --error="${WORKDIR}/DuneEngquistOsher2d_03-${SLURM_JOB_ID}.err" "${UGSHELL}" "${ARGS[@]}" -pre_smoothLim 2 -post_smoothLim 2 -riemman 3 &
PID6=$!

STATUS=0

wait "${PID1}" || { echo "ERROR: 8 MPI rank Rusanov failed."; STATUS=1; }
wait "${PID2}" || { echo "ERROR: 16 MPI rank Rusanov failed."; STATUS=1; }
wait "${PID3}" || { echo "ERROR: 24 MPI rank Rusanov failed."; STATUS=1; }
wait "${PID4}" || { echo "ERROR: 8 MPI rank Engquist-Osher failed."; STATUS=1; }
wait "${PID5}" || { echo "ERROR: 16 MPI rank Engquist-Osher failed."; STATUS=1; }
wait "${PID6}" || { echo "ERROR: 24 MPI rank Engquist-Osher failed."; STATUS=1; }



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
