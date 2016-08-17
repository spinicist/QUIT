#!/bin/bash -eu

# Tobias Wood 2015
# Simple test script for elliptical parameter mapping

source ./test_common.sh
SILENCE_TESTS="0"

DATADIR="es"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

# First, create input data
SIZE="32 32 35"
$QUITDIR/qinewimage PD.nii --size "$SIZE" -f 1.0
$QUITDIR/qinewimage T1.nii --size "$SIZE" -g "0 0.5 2.0"
$QUITDIR/qinewimage T2.nii --size "$SIZE" -g "1 0.01 0.25"
$QUITDIR/qinewimage f0.nii --size "$SIZE" -g "2 -150.0 150.0"
$QUITDIR/qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_TR="0.005"
SSFP_FLIP="15"
SSFP_PC4="0 90 180 270"
SSFP_PC6="0 60 120 180 240 300"

function run_tests() {
PREFIX="$1"
SSFP_PC="$2"
PREC="$3"
ARGS="$4"
NOISE="0.002"
run_test "CREATE_SIGNALS" $QUITDIR/qisignal --1 -n -x -v --noise $NOISE << END_SIG
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
$SSFP_TR" > ${PREFIX}es_in.txt

run_test "${PREFIX}"     $QUITDIR/qiesmap -n -v -bB1.nii ${PREFIX}${SSFP_FILE}  -o${PREFIX} $ARGS < ${PREFIX}es_in.txt
compare_test "${PREFIX}T1" T1.nii ${PREFIX}ES_T1.nii $NOISE $PREC
compare_test "${PREFIX}T2" T2.nii ${PREFIX}ES_T2.nii $NOISE $PREC
compare_test "${PREFIX}PD" PD.nii ${PREFIX}ES_M.nii $NOISE $PREC
compare_test "${PREFIX}f0" f0.nii ${PREFIX}ES_f0.nii $NOISE $PREC
}

run_tests "HYPER6" "$SSFP_PC6" "50" "-T1"

cd ..
SILENCE_TESTS="0"
