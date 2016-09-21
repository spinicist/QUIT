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
SIZE="16 16 101"
$QUITDIR/qinewimage PD.nii --size "$SIZE" -f "1.0"
$QUITDIR/qinewimage T1.nii --size "$SIZE" -g "0 1.0 2.0"
$QUITDIR/qinewimage T2.nii --size "$SIZE" -g "1 0.025 0.25"
$QUITDIR/qinewimage f0.nii --size "$SIZE" -g "2 -200.0 200.0"
$QUITDIR/qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
GS_FILE="gs.nii"
SSFP_TR="0.005"
SSFP_FLIP="45"
#SSFP_PC4="0 90 180 270"
SSFP_PC="180 240 300 0 60 120"

NOISE="0.002"
run_test "CREATE_SIGNALS" $QUITDIR/qisignal --model=1 -n -x -v << END_SIG
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
$GS_FILE
SSFP_GS
$SSFP_FLIP
$SSFP_TR
END
END_SIG

echo "$SSFP_FLIP
$SSFP_PC
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

run_test "GS"  $QUITDIR/qissfpbands ${SSFP_FILE} --ph_incs=6 --ph_order -o gs.nii
qicomplex -x gs.nii -M gs_mag.nii -P gs_ph.nii
qicomplex -x ssfp.nii -M ssfp_mag.nii -P ssfp_ph.nii
fslmaths gs_ph.nii -div 3.141592 -div $SSFP_TR gs_f0.nii
run_tests "HYPER" "50" "--algo=h --ph_incs=6 --ph_order"
run_tests "CONST" "50" "--algo=c --ph_incs=6 --ph_order"

cd ..
SILENCE_TESTS="0"
