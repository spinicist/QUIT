#!/bin/bash -e

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="mcdespot"
rm -rf $DATADIR
mkdir -p $DATADIR
cd $DATADIR

DIMS="5 5 4"

$QUITDIR/qnewimage -d "$DIMS" -f "1.0" PD.nii
$QUITDIR/qnewimage -d "$DIMS" -f "0.465" T1_a.nii
$QUITDIR/qnewimage -d "$DIMS" -f "0.026" T2_a.nii
$QUITDIR/qnewimage -d "$DIMS" -f "1.070" T1_b.nii
$QUITDIR/qnewimage -d "$DIMS" -f "0.117" T2_b.nii
$QUITDIR/qnewimage -d "$DIMS" -f "0.18" tau_a.nii
$QUITDIR/qnewimage -d "$DIMS" -g "0 -25. 25." f0.nii
$QUITDIR/qnewimage -d "$DIMS" -g "1 0.75 1.25" B1.nii
$QUITDIR/qnewimage -d "$DIMS" -g "2 0.1 0.25" f_a.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="3 4 5 6 7 9 13 18
0.0065"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 16 21 27 33 40 51 68"
SSFP_TR="0.005"
SSFP_PAR_180_0="$SSFP_FLIP
180 0
$SSFP_TR"
SSFP_PAR_90_270="$SSFP_FLIP
90 270
$SSFP_TR"

run_test "CREATE_SIGNALS" $QUITDIR/qsignal --2 -n << END_MCSIG
PD.nii
T1_a.nii
T2_a.nii
T1_b.nii
T2_b.nii
tau_a.nii
f_a.nii
f0.nii
B1.nii
SPGR
$SPGR_PAR
$SPGR_FILE
SSFP
$SSFP_PAR_180_0
180_0${SSFP_FILE}
SSFP
$SSFP_PAR_90_270
90_270${SSFP_FILE}
END
END_MCSIG

echo "SPGR #Signal type
# Test comment line
$SPGR_FILE # Filename
$SPGR_PAR
SSFP
180_0${SSFP_FILE}
$SSFP_PAR_180_0
END" > mcd.in

function run() {
PREFIX="$1"
OPTS="$2"
run_test $PREFIX $QUITDIR/qmcdespot $OPTS -o $PREFIX < mcd.in

echo "Tau:  " $( fslstats ${PREFIX}2C_tau_a.nii -m -s )
echo "T1_a: " $( fslstats ${PREFIX}2C_T1_a.nii -m -s )
echo "T2_a: " $( fslstats ${PREFIX}2C_T2_a.nii -m -s )
echo "T1_b: " $( fslstats ${PREFIX}2C_T1_b.nii -m -s )
echo "T2_b: " $( fslstats ${PREFIX}2C_T2_b.nii -m -s )

compare_test $PREFIX f_a.nii ${PREFIX}2C_f_a.nii 0.05

}

run "f0" "-v --2 -n -S1 -g -bB1.nii -ff0.nii "
#run "GAUSS" "-v --2 -n -S1 -g -bB1.nii"
#run "MEAN" "-v --2 -n -S1 -bB1.nii"


cd ..
SILENCE_TESTS="0"
