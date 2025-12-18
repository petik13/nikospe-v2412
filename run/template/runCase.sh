#!/bin/bash
topoSet
# renumberMesh -overwrite
decomposePar
foamJob -s -p renumberMesh -overwrite
nohup mpirun -np 32 PHIWaveCurSph2 -parallel -withFunctionObjects > log &
