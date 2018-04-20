#!/bin/bash
### USAGE: ./run-process.sh <n-points>

REPO=/lfs/l2/gerda/pertoldi/0nbb-sensitivity
a=0.00001
b=0.01

for i in `seq 1 $1`; do
    qsub -N znbb-sensitivity.$i $REPO/tools/qsub/process.qsub `echo "$a+($i-1)*($b-$a)/($1)" | bc -l`
done
