#!/bin/bash -e

# Tobias Wood 2015
# Simple test script for multiecho and similar programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="relax"
mkdir -p $DATADIR
cd $DATADIR

DIMS="11 11 11"
$QUITDIR/qnewimage -d "$DIMS" -f 1 PD.nii
$QUITDIR/qnewimage -d "$DIMS" -g "0 0.5 5" T1.nii
$QUITDIR/qnewimage -d "$DIMS" -g "1 0.05 0.5" T2.nii
$QUITDIR/qnewimage -d "$DIMS" -g "2 -25.0 25.0" f0.nii
$QUITDIR/qnewimage -d "$DIMS" -f 1 B1.nii
# Setup parameters
SPIN_FILE="multispinecho.nii"
SPIN_PAR="0.005
0.005
3"

# Create input for Single Component
MCSIG_INPUT="PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SPINECHO
$SPIN_PAR
$SPIN_FILE
END"
echo "$MCSIG_INPUT" > qsignal.in
run_test "CREATE_SIGNALS" $QUITDIR/qsignal --1 < qsignal.in

echo "$SPIN_PAR" > multiecho_spin.in

run_test "SPINECHO_LOGLIN" $QUITDIR/qmultiecho $SPIN_FILE -n -al -oL -v < multiecho_spin.in
run_test "SPINECHO_ARLO" $QUITDIR/qmultiecho $SPIN_FILE -n -aa -oA -v < multiecho_spin.in
run_test "SPINECHO_NONLIN" $QUITDIR/qmultiecho $SPIN_FILE -n -an -oN -v < multiecho_spin.in
compare_test "LOGLIN" T2.nii LME_T2.nii 0.01
compare_test "ARLO" T2.nii AME_T2.nii 0.01
compare_test "NONLIN" T2.nii NME_T2.nii 0.01


cd ..
SILENCE_TESTS="0"
