#!/bin/bash


ACML_DIR=/pica/h1/afshin/acml/acmllib/ifort64_fma4
ACML_LIB=${ACML_DIR}/lib/libacml.a
export LD_LIBRARY_PATH=${ACML_DIR}/lib:${LD_LIBRARY_PATH}


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
missed=../../../out/missed_time.py

set +x
module load python 
set -x
python $compare $ref $res>>$out

python $task_times $out > task_times.txt
python $missed $out >> $out

for f in $(ls *.dot 2>/dev/null)
do 
	dot $f -Tpng -o $f.png
done
