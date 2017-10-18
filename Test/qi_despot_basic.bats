# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "DESPOT1-Basic" {

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="3 3 20 20"
SPGR_TR="0.01"
SIZE="32,32,32"
NOISE="0.002"
qinewimage --size "$SIZE" -g "1 0.8 1.0" PD.nii
qinewimage --size "$SIZE" -g "0 0.5 1.5" T1.nii
qisignal --model=1 -v -n --noise=$NOISE << OUT
PD.nii
T1.nii



$SPGR_FILE
SPGR
$SPGR_FLIP
$SPGR_TR
END
OUT
qidespot1 -v -n $SPGR_FILE <<OUT
$SPGR_FLIP
$SPGR_TR
OUT
img_diff T1.nii D1_T1.nii $NOISE 50

}

@test "DESPOT2-Basic" {

# Setup parameters
FILE="ssfp.nii"
FLIP="15 60"
PINC="180"
TR="0.01"
SIZE="32,32,32"
NOISE="0.002"
qinewimage --size "$SIZE" -g "1 0.8 1.0" PD.nii
[ -e PD.nii ]
qinewimage --size "$SIZE" -g "0 0.5 1.5" T1.nii
[ -e T1.nii ]
qinewimage --size "$SIZE" -g "2 0.02 0.1" T2.nii
[ -e T2.nii ]
qisignal --model=1 -v -n --noise=$NOISE << OUT
PD.nii
T1.nii
T2.nii


$FILE
SSFP
$FLIP
$PINC
$TR
END
OUT
qidespot2 -v -n T1.nii $FILE <<OUT
$FLIP
$PINC
$TR
OUT
img_diff T2.nii D2_T2.nii $NOISE 50

}