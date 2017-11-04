#!/bin/sh
cd $1
rm *.dat
cp ../qmc .
mpirun -np 1 qmc  > qmc.out &
PID=$!
echo $PID
