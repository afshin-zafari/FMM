#!/bin/bash

p=$1
exp=${p}*/time*.txt
grep -ri "Mv\_far" $exp
grep -ri "Mv\_near" $exp
grep -ri "FMM" $exp
