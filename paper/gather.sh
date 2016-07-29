#!/bin/bash

base=../out/N1e5

P=050
Q=010
c1=$(ls $base/Exp*P${P}_Q${Q}*/trace*)
for f in $c1
do
   cp $f ./trace_files 
done
