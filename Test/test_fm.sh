#!/bin/bash

# Tobias Wood 2015
# Simple test script for FM, it runs slowly so test
# separately to other single-component programs.

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="0"

DATADIR="fm"
rm -rf $DATADIR
mkdir -p $DATADIR
cd $DATADIR

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
SSFP_FILE="ssfp.nii"
SSFP_FLIP="15 30 45"
SSFP_PC="0 180"
SSFP_TR="0.005"

echo "$MCSIG_INPUT" > signal.in
run_test "CREATE_REAL_SIGNALS" $QUITDIR/qisignal --1 -n -v << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SPGR
$SPGR_PAR
$SPGR_FILE
SSFP
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_FILE
END
END_SIG

echo "$SSFP_PAR" > despot2fm.in

run_test "FMLM" $QUITDIR/qidespot2fm T1.nii $SSFP_FILE -n -bB1.nii -v -o LM << END_FM
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END_FM

run_test "FMFINITE" $QUITDIR/qidespot2fm T1.nii $SSFP_FILE -n -bB1.nii -v -o F --finite << END_FM
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
0.001 #Dummy value for now
END_FM

run_test "FMSRC" $QUITDIR/qidespot2fm T1.nii $SSFP_FILE -n -as -bB1.nii -v -o SRC << END_FM

compare_test "FMLM" T2.nii LMFM_T2.nii 0.001
compare_test "FMFINITE" T2.nii FFM_T2.nii 0.001
compare_test "FMSRC" T2.nii SRCFM_T2.nii 0.01

cd ..
SILENCE_TESTS="0"
