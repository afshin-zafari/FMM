#!/bin/bash

#SBATCH -A p2009014
#SBATCH -o FMM-SCHED-DRAW-%j.out
#SBATCH -p devel
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J FMM-SCHED-DRAW


module load python/3.5.0

list=$(ls -R */*/*tra*_t_*.txt)

for f in $list
do
    if [ -e "${f}.pdf" ] ; then 
	echo "$f ----- Schedule drawn"
	continue;
    fi
    echo $f
    python drawsched.py $f
done


