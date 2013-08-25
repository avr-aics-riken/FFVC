#!/bin/sh

a=1
while [ $a -le 12 ]
do
export OMP_NUM_THREADS=$a
ffvc pm.tp
./change.sh $a
a=`expr $a + 1`
done
