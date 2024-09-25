#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -C cpu
#SBATCH --time=00:00:15
#SBATCH --qos=debug
#SBATCH --account=m4776

srun -n 2 -c 1 ./demo-submit.x
