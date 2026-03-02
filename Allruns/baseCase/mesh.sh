#!/bin/bash
set -e
surfaceFeatureExtract
blockMesh
decomposePar

foamJob -parallel -screen snappyHexMesh -overwrite

# Reconstruct the mesh into ./constant/polyMesh (root case)
reconstructParMesh -constant

# Optional: reconstruct fields too (if you ran/created any)
# reconstructPar -latestTime

# Delete processor dirs
rm -rf processor*
