#!/bin/bash -e

# Tobias Wood 2015
# Simple test script for FM, it runs slowly so test
# separately to other single-component programs.

# Tests whether programs run successfully on toy data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="fm"
rm -rf $DATADIR
mkdir -p $DATADIR
cd $DATADIR

# First, create input data
DIMS="10 10 21"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage PD.nii -d "$DIMS" -v "$VOXDIMS" -f 1.0
$QUITDIR/qinewimage T1.nii -d "$DIMS" -v "$VOXDIMS" -f 1.25
$QUITDIR/qinewimage T2.nii -d "$DIMS" -v "$VOXDIMS" -g "1 0.025 0.125"
$QUITDIR/qinewimage f0.nii -d "$DIMS" -v "$VOXDIMS" -g "2 -200.0 200.0"
$QUITDIR/qinewimage B1.nii -d "$DIMS" -v "$VOXDIMS" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_FLIP="10 15 20 30 40 50 60"
SSFP_NEG_FLIP="-15 -60"
SSFP_TR="0.006"
SSFP_Trf="0.002"

function run_tests() {
PREFIX="$1"
SSFP_PC="$2"

run_test "CREATE_SIGNALS" $QUITDIR/qisignal --1 -n -v -x --noise 0.002 << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SSFP
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
${PREFIX}x$SSFP_FILE
SSFPFinite
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_Trf
${PREFIX}xF${SSFP_FILE}
END
END_SIG

$QUITDIR/qicomplex -x ${PREFIX}x$SSFP_FILE -M ${PREFIX}$SSFP_FILE -P ${PREFIX}p$SSFP_FILE
$QUITDIR/qicomplex -x ${PREFIX}xF$SSFP_FILE -M ${PREFIX}F$SSFP_FILE -P ${PREFIX}pF$SSFP_FILE

echo "$SSFP_FLIP
$SSFP_PC
$SSFP_TR" > ${PREFIX}fm_in.txt

echo "$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_Trf" > ${PREFIX}fm_f_in.txt

run_test "FM"    $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}$SSFP_FILE    -o ${PREFIX}       < ${PREFIX}fm_in.txt
run_test "SRC"   $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}$SSFP_FILE    -as -o${PREFIX}s   < ${PREFIX}fm_in.txt
run_test "XFM"   $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}xF$SSFP_FILE   -ax -o${PREFIX}x   < ${PREFIX}fm_in.txt
run_test "FM_F"  $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}F${SSFP_FILE} --finite -o${PREFIX}F       < ${PREFIX}fm_f_in.txt
run_test "SRC_F" $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}F${SSFP_FILE} -as --finite -o${PREFIX}sF   < ${PREFIX}fm_f_in.txt
run_test "XFM_F" $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}xF${SSFP_FILE} --finite -ax -o${PREFIX}xF < ${PREFIX}fm_f_in.txt

compare_test "FM"    T2.nii ${PREFIX}FM_T2.nii   0.01
compare_test "SRC"   T2.nii ${PREFIX}sFM_T2.nii  0.01
compare_test "XFM"   T2.nii ${PREFIX}xFM_T2.nii  0.015
compare_test "FM_F"  T2.nii ${PREFIX}FFM_T2.nii  0.01
compare_test "SRC_F" T2.nii ${PREFIX}sFFM_T2.nii 0.01
compare_test "XFM_F" T2.nii ${PREFIX}xFFM_T2.nii 0.015
}

run_tests "2" "0 180"
run_tests "4" "0 90 180 270"

cd ..
SILENCE_TESTS="0"
