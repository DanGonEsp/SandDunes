#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=workq
#SBATCH -J SimpleAvalanche_2_Lev4
#SBATCH --output=/scratch/gonzald/SandDunes/SandDunes-output.txt
#SBATCH --error=/scratch/gonzald/SandDunes/SandDunes-error.txt
#SBATCH --mail-user=daniel.gonzalezesparza@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10105
#SBATCH -t 24:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1


#run the application 16:
srun -n 64 --ntasks-per-node=64 /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 2 -simCaseBnd 1 -timeMethod limex  -numProc 4 -numRefs 4 -numPreRefs 3 -numTimeSteps 100 -DT 10.0

#run the application 8:
srun -n 32 --ntasks-per-node=32 /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 2 -simCaseBnd 1 -timeMethod limex  -numProc 4 -numRefs 4 -numPreRefs 3 -numTimeSteps 100 -DT 10.0

#run the application 4:
srun -n 16 --ntasks-per-node=16 /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 2 -simCaseBnd 1 -timeMethod limex  -numProc 4 -numRefs 4 -numPreRefs 3 -numTimeSteps 100 -DT 10.0

#run the application 2:
srun -n 8 --ntasks-per-node=8 /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 2 -simCaseBnd 1 -timeMethod limex  -numProc 4 -numRefs 4 -numPreRefs 3 -numTimeSteps 100 -DT 10.0

#run the application 1:
srun -n 4 --ntasks-per-node=4 /project/k10105/gonzald/SandDunes/UG4/ug4/bin/ugshell -ex Avalanche.lua -dir_name /scratch/gonzald/SandDunes -dim 2 -simCaseBnd 1 -timeMethod limex  -numProc 4 -numRefs 4 -numPreRefs 3 -numTimeSteps 100 -DT 10.0







#End of File
