#!/bin/bash
rm -r 0
cp -r 0.orig 0
topoSet -dict system/topoSetDict
topoSet -dict system/topoSetDict_2
renumberMesh -overwrite
decomposePar
foamJob -s -p renumberMesh -overwrite
nohup mpirun -np 48 shipFlow -parallel -withFunctionObjects > log &
