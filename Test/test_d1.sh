#!/bin/bash -e

# Tobias Wood 2015
# Simple test scripts for QUIT programs

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="d1"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

DIMS="16 16 25"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 1.2 0.8" PD.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "0 0.5 2" T1.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "1 0.05 0.5" T2.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 -25.0 25.0" f0.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 0.75 1.25" B1.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="5 10 15 #Test comment
#Test line below has trailing whitespace
0.01 "
MPRAGE_FILE="mprage.nii"
MPRAGE_PAR="5
0.01
11
0.5
0"
HIFI_PAR="$SPGR_PAR
$MPRAGE_PAR"
AFI_FILE="afi.nii"
AFI_PAR="55.
0.02 0.1"
# Create input for Single Component
SIG_INPUT="PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SPGR_FILE
SPGR
$SPGR_PAR
$MPRAGE_FILE
MPRAGE
$MPRAGE_PAR
$AFI_FILE
AFI
$AFI_PAR
END"
echo "$SIG_INPUT" > signal.in
run_test "CREATE_REAL_SIGNALS" $QUITDIR/qisignal --1 -v -n < signal.in

echo "$SPGR_PAR" > despot1.in
echo "$HIFI_PAR" > despot1hifi.in

run_test "AFI" $QUITDIR/qiafi $AFI_FILE
compare_test "AFI_B1" B1.nii AFI_B1.nii 0.01
run_test "DESPOT1" $QUITDIR/qidespot1 $SPGR_FILE -n -bB1.nii < despot1.in
compare_test "DESPOT1" T1.nii D1_T1.nii 0.01
run_test "DESPOT1LM" $QUITDIR/qidespot1 $SPGR_FILE -n -an -oN -bB1.nii < despot1.in
compare_test "DESPOT1LM" T1.nii ND1_T1.nii 0.01
run_test "DESPOT1HIFI" $QUITDIR/qidespot1hifi $SPGR_FILE $MPRAGE_FILE -M -n -T1 < despot1hifi.in
compare_test "DESPOT1HIFI" T1.nii HIFI_T1.nii 0.01

cd ..
SILENCE_TESTS="0"
