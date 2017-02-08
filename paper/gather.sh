#!/bin/bash

base=../out4/N1e5

P=400
Q=100
c1=$(ls $base/Exp*P${P}_Q${Q}*/trace*)
for f in $c1
do
   cp $f ./trace_files 
done
