#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -C cpu
#SBATCH --time=00:00:15
#SBATCH --qos=debug
#SBATCH --account=m4776

srun -n 4 -c 1 ./hello-mpi.x
