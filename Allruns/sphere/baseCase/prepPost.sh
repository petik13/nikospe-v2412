#!/bin/bash
# Run after each case finishes: print the mean drift loads and motion RAOs,
# and append one row to the sweep file for collection across runs.

python3 meanLoads.py --append ../sweep_$(basename "$PWD" | sed 's/.*_head/head/').dat | tee summary.txt

reconstructPar > /dev/null 2>&1
rm -rf processor*
