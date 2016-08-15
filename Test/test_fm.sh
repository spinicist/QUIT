#!/bin/bash -eu

# Tobias Wood 2015
# Simple test script for FM, it runs slowly so test
# separately to other single-component programs.

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="fm"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

# First, create input data
SIZE="16 16 17"
$QUITDIR/qinewimage PD.nii --size "$SIZE" -f 1.0
$QUITDIR/qinewimage T1.nii --size "$SIZE" -f 1.0
$QUITDIR/qinewimage T2.nii --size "$SIZE" -g "1 0.025 0.125"
$QUITDIR/qinewimage f0.nii --size "$SIZE" -g "2 -150.0 150.0"
$QUITDIR/qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_TR="0.005"
SSFP_Trf="0.001"

function run_tests() {
PREFIX="$1"
SSFP_FLIP="$2"
SSFP_PC="$3"
PREC="$4"
ARGS="$5"
NOISE="0.002"
run_test "CREATE_SIGNALS" $QUITDIR/qisignal --1 -n -v --noise $NOISE << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
${PREFIX}$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
END_SIG

echo "$SSFP_FLIP
$SSFP_PC
$SSFP_TR" > ${PREFIX}fm_in.txt

run_test "${PREFIX}"     $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}${SSFP_FILE}  -o${PREFIX} $ARGS < ${PREFIX}fm_in.txt
compare_test "${PREFIX}"   T2.nii ${PREFIX}FM_T2.nii $NOISE $PREC

}

run_tests "BFGS2SYM" "12 65" "180 0" "50" ""
run_tests "BFGS2ASYM2" "12 65" "0 90 180 270" "50" "--asym"
run_tests "BFGS2ASYM6" "12 65" "0 90 180 270" "50" "--asym --off=6"
run_tests "CM2SYM" "12 65" "180 0" "50" "-ac"
run_tests "CM2ASYM" "12 65" "90 270" "50" "-ac --asym"

cd ..
SILENCE_TESTS="0"
