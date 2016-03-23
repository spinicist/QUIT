#!/bin/bash -e

# Tobias Wood 2015
# Simple test scripts for QUIT programs

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="2C"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

SIZE="20 5 3"
$QUITDIR/qinewimage --size "$SIZE" -f "1.0" PD.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.465" T1_m.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.026" T2_m.nii
$QUITDIR/qinewimage --size "$SIZE" -f "1.070" T1_ie.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.117" T2_ie.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.18" tau_m.nii
$QUITDIR/qinewimage --size "$SIZE" -g "0 0. 200." f0.nii
$QUITDIR/qinewimage --size "$SIZE" -g "1 0.75 1.25" B1.nii
$QUITDIR/qinewimage --size "$SIZE" -g "2 0.05 0.25" f_m.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="3 4 5 6 7 9 13 18"
SPGR_TR="0.0065"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 16 21 27 33 40 51 68 12 16 21 27 33 40 51 68"
SSFP_PC="180 180 180 180 180 180 180 180 0 0 0 0 0 0 0 0"
SSFP_TR="0.005"

INPUT="$SPGR_FILE
SPGR
$SPGR_FLIP
$SPGR_TR
$SSFP_FILE
SSFP_ECHO_FLEX
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
"

echo "$INPUT" > mcd_input.txt

run_test "CREATE_SIGNALS" $QUITDIR/qisignal --2 -v << END_SIG
PD.nii
T1_m.nii
T2_m.nii
T1_ie.nii
T2_ie.nii
tau_m.nii
f_m.nii
f0.nii
B1.nii
$INPUT
END_SIG

function run() {
PREFIX="$1"
OPTS="$2"
run_test $PREFIX $QUITDIR/qimcdespot $OPTS -M2 -bB1.nii -r -o $PREFIX -v -n < mcd_input.txt

compare_test $PREFIX f_m.nii ${PREFIX}2C_f_m.nii 0.05

#echo "       Mean     Std.     CoV"
#echo "T1_m:  " $( fslstats ${PREFIX}2C_T1_m.nii  -m -s | awk '{if(($1)>(0.)) {print $1, $2, $2/$1} else {print 0}}' )
#echo "T2_m:  " $( fslstats ${PREFIX}2C_T2_m.nii  -m -s | awk '{if(($1)>(0.)) {print $1, $2, $2/$1} else {print 0}}' )
#echo "T1_ie: " $( fslstats ${PREFIX}2C_T1_ie.nii -m -s | awk '{if(($1)>(0.)) {print $1, $2, $2/$1} else {print 0}}' )
#echo "T2_ie: " $( fslstats ${PREFIX}2C_T2_ie.nii -m -s | awk '{if(($1)>(0.)) {print $1, $2, $2/$1} else {print 0}}' )
#echo "MWF:   " $( fslstats ${PREFIX}2C_f_m.nii   -m -s | awk '{if(($1)>(0.)) {print $1, $2, $2/$1} else {print 0}}' )
#echo "Tau:   " $( fslstats ${PREFIX}2C_tau_m.nii -m -s | awk '{if(($1)>(0.)) {print $1, $2, $2/$1} else {print 0}}' )
}

run "BFGS" " -ff0.nii -ab"
run "SRC"  " -ff0.nii -aS"
run "GRC"  " -ff0.nii -aG"
cd ..
SILENCE_TESTS="0"
