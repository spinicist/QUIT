#!/bin/bash -eu

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

SIZE="32 32 32"
$QUITDIR/qinewimage --size "$SIZE" -g "1 0.75 1.25" PD.nii
$QUITDIR/qinewimage --size "$SIZE" -g "0 0.5 1.5" T1.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.05" T2.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.0" f0.nii
$QUITDIR/qinewimage --size "$SIZE" -g "2 0.75 1.25" B1.nii

# Setup parameters
SPGR2_FILE="spgr2.nii"
SPGR2_FLIP="3 3 20 20"
SPGR4_FILE="spgr4.nii"
SPGR4_FLIP="3 9 15 20"
SPGR_TR="0.01"
MPRAGE_FILE="mprage.nii"
MPRAGE_FLIP="5"
MPRAGE_SEGSIZE="64"
MPRAGE_KZERO="0"
MPRAGE_INV="0.45"
MPRAGE_DELAY="0.0"
AFI_FILE="afi.nii"
AFI_PAR="55.
0.02 0.1"
# Create input for Single Component
run_test "CREATE_REAL_SIGNALS" $QUITDIR/qisignal --1 -v -n --noise=0.002 <<END_SIGNALS
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SPGR2_FILE
SPGR
$SPGR2_FLIP
$SPGR_TR
$SPGR4_FILE
SPGR
$SPGR4_FLIP
$SPGR_TR
$MPRAGE_FILE
MPRAGE
$MPRAGE_FLIP
$SPGR_TR
$MPRAGE_SEGSIZE
$MPRAGE_KZERO
$MPRAGE_INV
$MPRAGE_DELAY
$AFI_FILE
AFI
$AFI_PAR
END
END_SIGNALS

echo "$SPGR2_FLIP
$SPGR_TR" > despot1_2.in
echo "$SPGR4_FLIP
$SPGR_TR" > despot1_4.in
echo "$SPGR2_FLIP
$SPGR_TR
$MPRAGE_FLIP
$SPGR_TR
$MPRAGE_SEGSIZE
$MPRAGE_KZERO
$MPRAGE_INV
$MPRAGE_DELAY" > despot1hifi.in

function d1_test () {
    TEST="$1"
    FILE="$2"
    ARGS="$3"
    INPUT="$4"
    run_test "DESPOT1_${TEST}" $QUITDIR/qidespot1 -v -n -bB1.nii $ARGS -o $TEST $FILE < $INPUT
    compare_test "DESPOT1_${TEST}" T1.nii ${TEST}D1_T1.nii 0.01
}

d1_test "LLS2"   $SPGR2_FILE "" despot1_2.in
d1_test "LLS4"   $SPGR4_FILE "" despot1_4.in
d1_test "LM"     $SPGR4_FILE "-an" despot1_4.in
#d1_test "LBFGSB" $SPGR4_FILE "-ab" despot1_4.in

run_test "HIFI" $QUITDIR/qidespot1hifi $SPGR2_FILE $MPRAGE_FILE -M -n -T1 -v < despot1hifi.in
compare_test "HIFIT1" T1.nii HIFI_T1.nii 0.01
compare_test "HIFIB1" B1.nii HIFI_B1.nii 0.05

run_test "AFI" $QUITDIR/qiafi $AFI_FILE -v
compare_test "AFI_B1" B1.nii AFI_B1.nii 0.05

cd ..
SILENCE_TESTS="0"
