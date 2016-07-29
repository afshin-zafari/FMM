#!/bin/bash

lstopo --force --of TXT topology.TXT
ACML_DIR=/pica/h1/afshin/acml/acmllib/ifort64_fma4
ACML_LIB=${ACML_DIR}/lib/libacml.a
export LD_LIBRARY_PATH=${ACML_DIR}/lib:${LD_LIBRARY_PATH}


expdir=Exp${j}_N${N}_P${S}_Q${Q}_${flags}_${T}cores
ex=$(ls -d Exp*_N${N}_P${S}_Q${Q}_${flags}_${T}cores)
curd=$(pwd)
if [ "z${ex}z" != "zz"   ] ; then 
    cp $ex/trace_f_t_S_w.txt $ex/trace_N${N}_P${S}_Q${Q}_${flags}_${T}cores
    return
    echo "----"
    cd $ex
else
    mkdir $expdir
    cd $expdir
fi




app=../../../bin/fmm

indir=../../../input
Tree="${indir}/tree_N${N}_L${L}_Q${Q}_S${S}.txt"
Ops="${indir}/operators_N${N}_L${L}_Q${Q}_S${S}.txt"
ref="${indir}/result_N${N}_L${L}_Q${Q}_S${S}.txt"
out=time_${flags}.txt

if [ ! -f $Tree ];
then
  echo "Not exist:    $Tree"
  cd $curd
  return
fi
$app $N $L $Q $S $T $Tree $Ops $flags >$out



res=$(ls ./result*)
compare=../../../out/check_results.py
task_times=../../../out/task_times.py
missed=../../../out/missed_time.py
dbrec=../../../out/dbrec.py

ttimes=task_times.txt

set +x
module load python 
set -x
python $compare $ref $res>>$out

python $task_times $out > $ttimes
python $missed $out >> $out
python $dbrec $N $L $S $Q $T $flags $out $ttimes >dbrec.txt

for f in $(ls *.dot 2>/dev/null)
do 
	dot $f -Tpng -o $f.png
done



cd $curd
