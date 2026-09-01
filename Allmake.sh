#!/bin/bash
#------------------------------------------------------------------------------
# Build all user libraries and solvers in this repository.
#
# Usage:
#   ./Allmake.sh              # wclean + wmake each target (full rebuild)
#   ./Allmake.sh --no-clean   # wmake only (incremental, much faster)
#
# Libraries are built before the solvers that link against them.
# A failing target does not stop the run; a summary is printed at the end
# and the script exits non-zero if anything failed.
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1                              # Run from this directory

DO_CLEAN=1
case "$1" in
    --no-clean|-n) DO_CLEAN=0 ;;
    "") ;;
    *) echo "Unknown option: $1"; echo "Usage: $0 [--no-clean]"; exit 1 ;;
esac

targets=(
    # --- Libraries (built first; the solvers below link against these) ---
    src/myFiniteVolume                                                    # libmyMovingWallSlip
    src/myFvMotionSolver                                                  # libmyFvMotionSolvers
    src/finiteVolume/fields/fvPatchFields/derived/myAdvective             # libmyAdvective
    src/finiteVolume/fields/fvPatchFields/derived/rigidTurgutBody         # librigidTurgutBody
    src/finiteVolume/fields/fvPatchFields/derived/linearizedRigidBody     # liblinearizedRigidBody
    src/finiteVolume/fields/fvPatchFields/derived/linBodyMotion           # liblinBodyMotion
    src/finiteVolume/fields/fvPatchFields/derived/waveCurrentPotential3D  # libwaveCurrentPotential3D
    src/finiteVolume/fields/fvPatchFields/derived/potForwardSpeedBC       # libpotForwardSpeedBC
    src/functionObjects/forces/myFunctionObject                           # libmyFunctionObject
    src/functionObjects/forces/myMeanForce                                # libmyMeanForce
    src/functionObjects/forces/meanWaveLoads                              # libMeanWaveLoads
    src/functionObjects/forces/middleFieldForm                            # libmiddleFieldForm

    # --- Solvers ---
    applications/solvers/PHIWaveCurSph2                                   # PHIWaveCurSph2
    applications/solvers/shipFlow                                         # shipFlow
    applications/solvers/manFlow                                          # manFlow
    applications/solvers/incompressible/myPimpleFoam                      # myPimpleFoam
    applications/solvers/incompressible/mySimpleFoam                      # mySimpleFoam
)

#------------------------------------------------------------------------------
# Deliberately NOT built:
#
#   src/finiteVolume/fields/fvPatchFields/derived/OLDlinBodyMotion
#       Archived copy. Its Make/files builds the SAME target as linBodyMotion
#       ($(FOAM_USER_LIBBIN)/liblinBodyMotion), so building both would mean one
#       silently overwrites the other depending on build order.
#------------------------------------------------------------------------------

failed=()
built=0
root="$PWD"

for dir in "${targets[@]}"; do
    if [ ! -f "$root/$dir/Make/files" ]; then
        echo "==> SKIP (no Make/files): $dir"
        failed+=("$dir (missing Make/files)")
        continue
    fi

    echo "=============================================================="
    echo "==> $dir"
    echo "=============================================================="

    cd "$root/$dir" || { failed+=("$dir (cd failed)"); continue; }

    [ "$DO_CLEAN" -eq 1 ] && wclean

    if wmake; then
        built=$((built + 1))
    else
        failed+=("$dir")
    fi

    cd "$root" || exit 1
done

#------------------------------------------------------------------------------
echo
echo "=============================================================="
echo "Built OK : $built / ${#targets[@]}"
if [ ${#failed[@]} -gt 0 ]; then
    echo "FAILED   : ${#failed[@]}"
    for f in "${failed[@]}"; do echo "   - $f"; done
    echo "=============================================================="
    exit 1
fi
echo "All targets built successfully."
echo "=============================================================="
