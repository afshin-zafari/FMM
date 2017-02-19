#!/bin/bash

for f in $(ls Exp*/dbrec.txt)
do
  cat $f
  cat $f >> alldbrec.txt
done