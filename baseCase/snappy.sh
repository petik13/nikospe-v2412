#!/bin/bash
set -e

# decomposePar

foamJob -parallel -screen snappyHexMesh -overwrite

# Reconstruct the mesh into ./constant/polyMesh (root case)
reconstructParMesh -constant

# Delete processor dirs
rm -rf processor*
