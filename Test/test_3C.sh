#!/bin/bash -e

# Tobias Wood 2015
# Simple test scripts for QUIT programs

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="3C"
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
$QUITDIR/qinewimage --size "$SIZE" -f 4.0 T1_csf.nii
$QUITDIR/qinewimage --size "$SIZE" -f 2.5 T2_csf.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.18" tau_m.nii
$QUITDIR/qinewimage --size "$SIZE" -g "0 -5. 5." f0.nii
$QUITDIR/qinewimage --size "$SIZE" -g "1 0.75 1.25" B1.nii
$QUITDIR/qinewimage --size "$SIZE" -g "2 0.05 0.25" f_m.nii
$QUITDIR/qinewimage --size "$SIZE" -f 0.05 f_csf.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="3 4 5 6 7 9 13 18"
SPGR_TR="0.0065"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 16 21 27 33 40 51 68 "
SSFP_PC="180 0"
SSFP_TR="0.005"

INPUT="$SPGR_FILE
SPGR
$SPGR_FLIP
$SPGR_TR
$SSFP_FILE
SSFP
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
"

echo "$INPUT" > mcd_input.txt

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
$INPUT
END_MCSIG

function run() {
PREFIX="$1"
OPTS="$2"
run_test $PREFIX $QUITDIR/qimcdespot $OPTS -n -v -M3 -bB1.nii -r -o $PREFIX < mcd_input.txt
compare_test $PREFIX f_m.nii ${PREFIX}3C_f_m.nii 0.05
}

run "GRC" "-ff0.nii -aG"

#echo "       Mean     Std.     CoV"
#echo "T1_m:  " $( fslstats ${PREFIX}${COMP}C_T1_m.nii  -m -s | awk '{print $1, $2, $2/$1}' )
#echo "T2_m:  " $( fslstats ${PREFIX}${COMP}C_T2_m.nii  -m -s | awk '{print $1, $2, $2/$1}' )
#echo "T1_ie: " $( fslstats ${PREFIX}${COMP}C_T1_ie.nii -m -s | awk '{print $1, $2, $2/$1}' )
#echo "T2_ie: " $( fslstats ${PREFIX}${COMP}C_T2_ie.nii -m -s | awk '{print $1, $2, $2/$1}' )
#echo "MWF:   " $( fslstats ${PREFIX}${COMP}C_f_m.nii   -m -s | awk '{print $1, $2, $2/$1}' )
#echo "Tau:   " $( fslstats ${PREFIX}2C_tau_m.nii -m -s | awk '{print $1, $2, $2/$1}' )


cd ..
SILENCE_TESTS="0"
