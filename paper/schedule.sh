#!/bin/bash

tf=$(ls ./trace_graphs/*/trace*SF_*)
tn=$(ls ./trace_graphs/*/trace*SN_*)
tm=$(ls ./trace_graphs/*/trace*SNF_*)

for f in $tn
do
   p="#ff0000"
   python drawsched.py $f $p
done

for f in $tf
do
   p="#0000ff #808080 #000000 #800080 #008000"
   python drawsched.py $f $p
done

for f in $tm
do
   p="#ff0000 #0000ff #808080 #000000 #800080 #008000"
   python drawsched.py $f $p
done
