#!/bin/bash
cd applications/solvers/PHIWaveCurSph2
wclean
wmake
cd -
cd src/functionObjects/forces/meanWaveLoads/
wclean
wmake
cd -
cd src/finiteVolume/fields/fvPatchFields/derived/waveCurrentPotential3D/
wclean
wmake
cd -
cd src/functionObjects/forces/myMeanForce/
wclean
wmake
cd -
cd src/functionObjects/forces/myFunctionObject/
wclean
wmake
cd -
cd src/finiteVolume/fields/fvPatchFields/derived/rigidTurgutBody/
wclean
wmake
cd -