#!/bin/bash -e

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
$QUITDIR/qinewimage T1.nii --size "$SIZE" -f 1.25
$QUITDIR/qinewimage T2.nii --size "$SIZE" -g "1 0.025 0.125"
$QUITDIR/qinewimage f0.nii --size "$SIZE" -g "2 -150.0 150.0"
$QUITDIR/qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_TR="0.01"
SSFP_Trf="0.001"

function run_tests() {
PREFIX="$1"
SSFP_FLIP="$2"
SSFP_PC="$3"

run_test "CREATE_SIGNALS" $QUITDIR/qisignal --1 -n -v -x --noise 0.002 << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
${PREFIX}x$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
${PREFIX}xf${SSFP_FILE}
SSFP_FINITE
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_Trf
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

run_test "${PREFIX}LM"       $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}${SSFP_FILE}  -o${PREFIX}   -al          < ${PREFIX}fm_in.txt
#run_test "${PREFIX}LM_F"     $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}f${SSFP_FILE} -o${PREFIX}f  -al --finite < ${PREFIX}fm_f_in.txt
run_test "${PREFIX}BFGS"     $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}${SSFP_FILE}  -o${PREFIX}b  -ab          < ${PREFIX}fm_in.txt
#run_test "${PREFIX}BFGS_F"   $QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ${PREFIX}f${SSFP_FILE} -o${PREFIX}bf -ab --finite < ${PREFIX}fm_f_in.txt

compare_test "LM"     T2.nii ${PREFIX}FM_T2.nii  0.01
#compare_test "LM_F"   T2.nii ${PREFIX}fFM_T2.nii 0.01
compare_test "BFGS"   T2.nii ${PREFIX}bFM_T2.nii  0.01
#compare_test "BFGS_F" T2.nii ${PREFIX}bfFM_T2.nii 0.01
}

run_tests "2" "5 25 45 65 5 25 45 65" "180 180 180 180 0 0 0 0"
#run_tests "4" "0 90 180 270"

cd ..
SILENCE_TESTS="0"
