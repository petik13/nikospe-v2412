#!/bin/bash
cd applications/solvers/PHIWaveCurSph2
wclean
wmake
cd -
cd src/functionObjects/forces/meanWaveLoads/
wclean
wmake
cd -
cd src/finiteVolume/fields/fvPatchFields/derived/waveCurMqLeastAB3DPotUPFD5/
wclean
wmake
cd -
cd src/functionObjects/forces/myMeanForce/
wclean
wmake
cd -