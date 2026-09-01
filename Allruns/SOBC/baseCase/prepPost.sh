#!/bin/bash
# Run after each case finishes: print the mean drift loads and motion RAOs,
# append one row to the per-condition sweep file, and collect the drift loads
# for the whole sweep into a single csv.
#
# The .dat name carries BOTH speed and heading (lam6.0767_U0.275_head30 ->
# sweep_U0.275_head30.dat).  Rows there are keyed on lambda/L, so a file that
# mixed speeds or headings would have every case overwrite the last.

cond=$(basename "$PWD" | sed 's/^lam[^_]*_//')

python3 meanLoads.py \
    --append "../sweep_${cond}.dat" \
    --csv    "../results.csv" | tee summary.txt

reconstructPar > /dev/null 2>&1
rm -rf processor*
