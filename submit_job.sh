#!/bin/bash
#SBATCH --job-name=1_singleion
#SBATCH --account=PAA0027
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --time=60:00:00
#SBATCH --output=1_jobinfo_%j.out

cd /users/PAA0027/yuanhaozhang/yz96/1_Tg/asiftest


module purge 
module load gcc/12.3.0
module load mvapich/4.1
module load python/3.12
python simulator.py
