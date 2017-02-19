#!/bin/bash

#SBATCH -A p2009014
#SBATCH -o FMM-N1e5-%j.out
#SBATCH -p node
#SBATCH -t 00:50:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J FMM-N1e5


set -x


N=100000
j=${SLURM_JOBID}

for F in MF MFN MN
do
    flags=ftwS$F
    for S in 050 300 400 
    do
	for Q in 010 020 050 100
	do
	    if [ "$S" -eq "400" ] ; then L=06; fi
	    if [ "$S" -eq "200" ] ; then L=07; fi
	    if [ "$S" -eq "300" ] ; then L=07; fi
	    if [ "$S" -eq "100" ] ; then L=08; fi
	    if [ "$S" -eq "050" ] ; then L=09; fi
	    for T in 1 4 8 16
	    do
	       OMP_NUM_THREADS=$T
	       export OMP_NUM_THREADS=$T
		source ../main.sh
	    done
	    T=1;
	    flags=fswS$F
	    source ../main.sh
	    flags=ftwS$F
	done
    done
done
