#!/bin/bash -eu

# Tobias Wood 2015
# Simple test script for elliptical parameter mapping

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="es"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

# First, create input data
SIZE="16 16 17"
$QUITDIR/qinewimage PD.nii --size "$SIZE" -g "2 0.9 1.1"
$QUITDIR/qinewimage T1.nii --size "$SIZE" -g "0 0.5 2.0"
$QUITDIR/qinewimage T2.nii --size "$SIZE" -g "1 0.01 0.25"
$QUITDIR/qinewimage f0.nii --size "$SIZE" -g "2 -150.0 150.0"
$QUITDIR/qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_TR="0.005"
SSFP_FLIP="20"
#SSFP_PC4="0 90 180 270"
SSFP_PC="0 30 60 90 120 150 180 210 240 270 300 330"

NOISE="0.002"
run_test "CREATE_SIGNALS" $QUITDIR/qisignal --1 -n -x -v --noise $NOISE << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
END_SIG

echo "$SSFP_FLIP
$SSFP_TR" > es_in.txt

function run_tests() {
PREFIX="$1"
PREC="$2"
ARGS="$3"

run_test "${PREFIX}"     $QUITDIR/qiesmap -n -v -bB1.nii ${SSFP_FILE}  -o${PREFIX} $ARGS < es_in.txt

compare_test "${PREFIX}T1" T1.nii ${PREFIX}ES_T1.nii $NOISE $PREC
compare_test "${PREFIX}T2" T2.nii ${PREFIX}ES_T2.nii $NOISE $PREC
compare_test "${PREFIX}PD" PD.nii ${PREFIX}ES_M.nii $NOISE $PREC
compare_test "${PREFIX}f0" f0.nii ${PREFIX}ES_f0.nii $NOISE $PREC

}

run_tests "CONST" "50" "--algo=c --ph_incs=12 --ph_order"
run_test "GS"  $QUITDIR/qissfpbands ${SSFP_FILE} --ph_incs=6 --ph_order --magnitude -o gs.nii
run_tests "HYPER" "50" "--algo=h --ph_incs=12 --ph_order"

cd ..
SILENCE_TESTS="0"
