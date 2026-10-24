
#!/bin/bash -l

#SBATCH --job-name=GMG-Benchmark
#SBATCH --partition=workq
#SBATCH --account=k10105
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=01:00:00
#SBATCH --array=0-7%8
#SBATCH --exclusive
#SBATCH --output=/scratch/gonzald/SandDunes/GMG-%A-%a.out
#SBATCH --error=/scratch/gonzald/SandDunes/GMG-%A-%a.err
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=FAIL,END

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

UGSHELL=/project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell
WORKDIR=/scratch/gonzald/SandDunes
APPDIR=/project/k10105/gonzald/SandDunes/Files

cd "${APPDIR}" || exit 1

BASE=(5 5 5 5 4 4 4 4)
PRE=(1 1 2 2 1 1 2 2)
POST=(0 1 0 1 0 1 0 1)

ID=${SLURM_ARRAY_TASK_ID}

CASE_DIR="${WORKDIR}/GMG_${SLURM_ARRAY_JOB_ID}_${ID}"
mkdir -p "${CASE_DIR}" || exit 1

echo "=============================================="
echo "GMG benchmark ${ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Base level: ${BASE[$ID]}"
echo "Pre-smoothing: ${PRE[$ID]}"
echo "Post-smoothing: ${POST[$ID]}"
echo "=============================================="

ARGS=(
    -ex Multiphase.lua
    -dir_name "${CASE_DIR}"
    -dim 2
    -timeMethod limex
    -numRefs 7
    -numPreRefs 3
    -numTimeSteps 100
    -DT 3600.0
    -gmgBaseLevel "${BASE[$ID]}"
    -pre_smoothImp "${PRE[$ID]}"
    -post_smoothImp "${POST[$ID]}"
    -max_newton_steps_steady_state 20
    -precondImp gmg
    -boolLoadCheckPoint false
    -boolSaveCheckPoint false
    -steadyStateOnly true
)

srun --exclusive --exact \
    --ntasks=8 \
    --ntasks-per-node=8 \
    --cpus-per-task=1 \
    --hint=nomultithread \
    --cpu-bind=cores \
    --mem=0 \
    "${UGSHELL}" "${ARGS[@]}"

STATUS=$?

if [ "${STATUS}" -ne 0 ]; then
    echo "ERROR: GMG benchmark ${ID} failed."
    exit "${STATUS}"
fi

echo "=============================================="
echo "GMG benchmark ${ID} finished."
echo "Job ID: ${SLURM_JOB_ID}"
echo "=============================================="
