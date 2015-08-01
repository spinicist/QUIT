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
DIMS="5 10 31"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage PD.nii -d "$DIMS" -v "$VOXDIMS" -f 1.0
$QUITDIR/qinewimage T1.nii -d "$DIMS" -v "$VOXDIMS" -f 1.25
$QUITDIR/qinewimage T2.nii -d "$DIMS" -v "$VOXDIMS" -g "1 0.025 0.125"
$QUITDIR/qinewimage f0.nii -d "$DIMS" -v "$VOXDIMS" -g "2 -150.0 150.0"
$QUITDIR/qinewimage B1.nii -d "$DIMS" -v "$VOXDIMS" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_FLIP="10 15 20 30 40 50 60"
SSFP_TR="0.01"
SSFP_Trf="0.001"

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
${PREFIX}xf${SSFP_FILE}
END
END_SIG

$QUITDIR/qicomplex -x ${PREFIX}x$SSFP_FILE -M ${PREFIX}$SSFP_FILE -P ${PREFIX}p$SSFP_FILE
$QUITDIR/qicomplex -x ${PREFIX}xf$SSFP_FILE -M ${PREFIX}f$SSFP_FILE -P ${PREFIX}pf$SSFP_FILE

echo "$SSFP_FLIP
$SSFP_PC
$SSFP_TR" > ${PREFIX}fm_in.txt

echo "$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_Trf" > ${PREFIX}fm_f_in.txt

run_test "LM"       $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii -T1 ${PREFIX}${SSFP_FILE}  -o${PREFIX}   -al          < ${PREFIX}fm_in.txt
run_test "LM_F"     $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii -T1 ${PREFIX}f${SSFP_FILE} -o${PREFIX}f  -al --finite < ${PREFIX}fm_f_in.txt
run_test "BFGS"     $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii -T1 ${PREFIX}${SSFP_FILE}  -o${PREFIX}b  -ab          < ${PREFIX}fm_in.txt
run_test "BFGS_F"   $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii -T1 ${PREFIX}f${SSFP_FILE} -o${PREFIX}bf -ab --finite < ${PREFIX}fm_f_in.txt

compare_test "LM"     T2.nii ${PREFIX}FM_T2.nii  0.01
compare_test "LM_F"   T2.nii ${PREFIX}fFM_T2.nii 0.01
compare_test "BFGS"   T2.nii ${PREFIX}bFM_T2.nii  0.01
compare_test "BFGS_F" T2.nii ${PREFIX}bfFM_T2.nii 0.01
}

run_tests "2" "0 180"
run_tests "4" "0 90 180 270"

cd ..
SILENCE_TESTS="0"
