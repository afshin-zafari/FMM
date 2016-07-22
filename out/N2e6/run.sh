#!/bin/bash

#SBATCH -A p2009014
#SBATCH -o FMM-N1e6-%j.out
#SBATCH -p devel
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J FMM-N1e6


set -x


j=${SLURM_JOBID}
flags=ftaS
T=16
N=1000000
L=12
Q=100
S=100


source ../main.sh
