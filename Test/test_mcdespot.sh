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

DIMS="11 11 4"

$QUITDIR/qinewimage -d "$DIMS" -f "1.0" PD.nii
$QUITDIR/qinewimage -d "$DIMS" -f "0.465" T1_m.nii
$QUITDIR/qinewimage -d "$DIMS" -f "0.026" T2_m.nii
$QUITDIR/qinewimage -d "$DIMS" -f "1.070" T1_ie.nii
$QUITDIR/qinewimage -d "$DIMS" -f "0.117" T2_ie.nii
$QUITDIR/qinewimage -d "$DIMS" -f "0.18" tau_m.nii
$QUITDIR/qinewimage -d "$DIMS" -g "0 -5. 5." f0.nii
$QUITDIR/qinewimage -d "$DIMS" -g "1 0.75 1.25" B1.nii
$QUITDIR/qinewimage -d "$DIMS" -g "2 0.1 0.25" f_m.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="3 4 5 6 7 9 13 18
0.0065"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 16 21 27 33 40 51 68 "
SSFP_TR="0.005"
SSFP_PAR_180_0="$SSFP_FLIP
180 0
$SSFP_TR"
SSFP_PAR_90_270="$SSFP_FLIP
90 270
$SSFP_TR"

run_test "CREATE_SIGNALS" $QUITDIR/qisignal --2 -n << END_MCSIG
PD.nii
T1_m.nii
T2_m.nii
T1_ie.nii
T2_ie.nii
tau_m.nii
f_m.nii
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
run_test $PREFIX $QUITDIR/qimcdespot $OPTS -2 -n -bB1.nii -r -o $PREFIX -v < mcd.in

echo "       Mean     Std.     CoV"
echo "T1_m:  " $( fslstats ${PREFIX}2C_T1_m.nii -m -s | awk '{print $1, $2, $2/$1}' )
echo "T2_m:  " $( fslstats ${PREFIX}2C_T2_m.nii -m -s | awk '{print $1, $2, $2/$1}' )
echo "T1_ie: " $( fslstats ${PREFIX}2C_T1_ie.nii -m -s | awk '{print $1, $2, $2/$1}' )
echo "T2_ie: " $( fslstats ${PREFIX}2C_T2_ie.nii -m -s | awk '{print $1, $2, $2/$1}' )
echo "MWF:   " $( fslstats ${PREFIX}2C_f_m.nii -m -s | awk '{print $1, $2, $2/$1}' )
echo "Tau:   " $( fslstats ${PREFIX}2C_tau_m.nii -m -s | awk '{print $1, $2, $2/$1}' )
compare_test $PREFIX f_m.nii ${PREFIX}2C_f_m.nii 0.05
}

run "GAUSSSCALE1f0" " -S1 -ff0.nii"
#run "SCALE1f0"      " -S1 -g0 -ff0.nii"
#run "NOSCALEf0"     " -SNONE -g0 -ff0.nii"
#run "MEANSCALEf0"   " -g0 -ff0.nii"
#run "MEANSCALE"     " -g0 "
#run "GAUSSSCALE1"   " -S1 -g0 "
#run "NOSCALE"       " -SNONE -g0 "

cd ..
SILENCE_TESTS="0"
