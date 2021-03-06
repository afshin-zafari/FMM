#!/bin/bash

#SBATCH -A p2009014
#SBATCH -o FMM-SCHED-DRAW-%j.out
#SBATCH -p devel
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J FMM-SCHED-DRAW


module load python/3.5.0

list=$(ls -R */*/*prog_dur*.txt)

for f in $list
do
    cat $f >> all_dur.txt
done


