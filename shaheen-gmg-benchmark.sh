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






export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

UGSHELL=/project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell
WORKDIR=/lustre2/project/k10105/gonzald/SandDunes/Files

cd "$WORKDIR" || exit 1

BASE=(5 5 5 5 4 4 4 4)
PRE=(1 1 2 2 1 1 2 2)
POST=(0 1 0 1 0 1 0 1)

ID=$SLURM_ARRAY_TASK_ID

echo "Test: $ID"
echo "Base: ${BASE[$ID]}"
echo "Pre: ${PRE[$ID]}"
echo "Post: ${POST[$ID]}"

srun --mpi=pmix "$UGSHELL" -ex Multiphase.lua \
    -DT 3600 \
    -numRefs 7 \
    -numPreRefs 3 \
    -gmgBaseLevel "${BASE[$ID]}" \
    -pre_smooth "${PRE[$ID]}" \
    -post_smooth "${POST[$ID]}" \
    -max_newton_steps_steady_state 20 \
    -precondImp gmg \
    -boolLoadCheckPoint false \
    -boolSaveCheckPoint false \
    -dir_name "GMG_Test_${SLURM_ARRAY_JOB_ID}_${ID}"
