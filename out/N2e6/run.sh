#!/bin/bash

#SBATCH -A p2009014
#SBATCH -o FMM-N2e6-%j.out
#SBATCH -p node
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J FMM-N2e6


set -x


j=${SLURM_JOBID}
flags=ftwSMNF
T=16
expdir=Exp${j}_${flags}_${T}cores
mkdir $expdir
cd $expdir

app=../../../bin/fmm
N=2000000
L=15
indir=../../../input
out=time_${flags}.txt

$app $N $L $flags $indir $T >$out


ref=$indir/results_2000000.txt
res=$(ls ./result*)
compare=../../../out/check_results.py
module load python 
python $compare $ref $res>>$out
