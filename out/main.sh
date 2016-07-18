#!/bin/bash


expdir=Exp${j}_${flags}_${T}cores
mkdir $expdir
cd $expdir

app=../../../bin/fmm

indir=../../../input
Tree="${indir}/tree_N${N}_L${L}_Q${Q}_S${S}.txt"
Ops="${indir}/operators_N${N}_L${L}_Q${Q}_S${S}.txt"
ref="${indir}/result_N${N}_L${L}_Q${Q}_S${S}.txt"
out=time_${flags}.txt

$app $N $L $Q $S $T $Tree $Ops $flags >$out



res=$(ls ./result*)
compare=../../../out/check_results.py
task_times=../../../out/task_times.py
set +x
module load python 
set -x
python $compare $ref $res>>$out

python $task_times $out > task_times.txt
