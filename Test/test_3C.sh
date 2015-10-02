#!/bin/bash -e

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="3C"
mkdir -p $DATADIR
cd $DATADIR
rm *

DIMS="5 5 3"

$QUITDIR/qinewimage -d "$DIMS" -f 1.0 PD.nii
$QUITDIR/qinewimage -d "$DIMS" -f 0.465 T1_m.nii
$QUITDIR/qinewimage -d "$DIMS" -f 0.026 T2_m.nii
$QUITDIR/qinewimage -d "$DIMS" -f 1.070 T1_ie.nii
$QUITDIR/qinewimage -d "$DIMS" -f 0.117 T2_ie.nii
$QUITDIR/qinewimage -d "$DIMS" -f 4.0 T1_csf.nii
$QUITDIR/qinewimage -d "$DIMS" -f 2.5 T2_csf.nii
$QUITDIR/qinewimage -d "$DIMS" -f "0.18" tau_m.nii
$QUITDIR/qinewimage -d "$DIMS" -g "0 -5. 5." f0.nii
$QUITDIR/qinewimage -d "$DIMS" -g "1 0.75 1.25" B1.nii
$QUITDIR/qinewimage -d "$DIMS" -g "2 0.05 0.25" f_m.nii
$QUITDIR/qinewimage -d "$DIMS" -f 0.05 f_csf.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 "
SPGR_TR="0.0065"
#SPGR_Trf="0.001"
#SPGR_TE="0.003"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="5 10 15 20 25 30 35 40 45 50 55 60 65 "
SSFP_PC="180 0"
SSFP_TR="0.005"
#SSFP_Trf="0.0025"

run_test "CREATE_SIGNALS" $QUITDIR/qisignal --3 -n -v << END_MCSIG
PD.nii
T1_m.nii
T2_m.nii
T1_ie.nii
T2_ie.nii
T1_csf.nii
T2_csf.nii
tau_m.nii
f_m.nii
f_csf.nii
f0.nii
B1.nii
SPGR
$SPGR_FLIP
$SPGR_TR
$SPGR_FILE
SSFP
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_FILE
END
END_MCSIG

function run() {
PREFIX="$1"
OPTS="$2"
run_test $PREFIX $QUITDIR/qimcdespot $OPTS -n -v -M3 -bB1.nii -r -o $PREFIX << END_INPUT
SPGR
$SPGR_FILE
$SPGR_FLIP
$SPGR_TR
SSFP
$SSFP_FILE
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
END_INPUT

compare_test $PREFIX f_m.nii ${PREFIX}3C_f_m.nii 0.05
}

run "GRC" "-ff0.nii -ab"

#echo "       Mean     Std.     CoV"
#echo "T1_m:  " $( fslstats ${PREFIX}${COMP}C_T1_m.nii  -m -s | awk '{print $1, $2, $2/$1}' )
#echo "T2_m:  " $( fslstats ${PREFIX}${COMP}C_T2_m.nii  -m -s | awk '{print $1, $2, $2/$1}' )
#echo "T1_ie: " $( fslstats ${PREFIX}${COMP}C_T1_ie.nii -m -s | awk '{print $1, $2, $2/$1}' )
#echo "T2_ie: " $( fslstats ${PREFIX}${COMP}C_T2_ie.nii -m -s | awk '{print $1, $2, $2/$1}' )
#echo "MWF:   " $( fslstats ${PREFIX}${COMP}C_f_m.nii   -m -s | awk '{print $1, $2, $2/$1}' )
#echo "Tau:   " $( fslstats ${PREFIX}2C_tau_m.nii -m -s | awk '{print $1, $2, $2/$1}' )


cd ..
SILENCE_TESTS="0"
