# Copyright Tobias Wood 2017
# Tests for mcDESPOT

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "2C mcDESPOT" {

SIZE="4,4,4"
qinewimage --size "$SIZE" -f "1.0" PD.nii
qinewimage --size "$SIZE" -f "0.465" T1_m.nii
qinewimage --size "$SIZE" -f "0.026" T2_m.nii
qinewimage --size "$SIZE" -f "1.070" T1_ie.nii
qinewimage --size "$SIZE" -f "0.117" T2_ie.nii
qinewimage --size "$SIZE" -f "0.18" tau_m.nii
qinewimage --size "$SIZE" -g "0 0. 200." f0.nii
qinewimage --size "$SIZE" -g "1 0.75 1.25" B1.nii
qinewimage --size "$SIZE" -g "2 0.05 0.25" f_m.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="3 4 5 6 7 9 13 18"
SPGR_TR="0.0065"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 16 21 27 33 40 51 68"
SSFP_PC="180 0"
SSFP_TR="0.005"

INPUT="$SPGR_FILE
SPGR
$SPGR_FLIP
$SPGR_TR
$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
"

echo "$INPUT" > mcd_input.txt
NOISE="0.002"
qisignal --model=2 -v --noise=$NOISE<< END_SIG
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

qimcdespot $OPTS -M2 -bB1.nii -ff0.nii -r -v -n < mcd_input.txt
qidiff --baseline=f_m.nii --input=2C_f_m.nii --noise=$NOISE --tolerance=250 --verbose

}

@test "3C mcDESPOT" {

SIZE="4,4,4"
qinewimage --size "$SIZE" -f "1.0" PD.nii
qinewimage --size "$SIZE" -f "0.465" T1_m.nii
qinewimage --size "$SIZE" -f "0.026" T2_m.nii
qinewimage --size "$SIZE" -f "1.070" T1_ie.nii
qinewimage --size "$SIZE" -f "0.117" T2_ie.nii
qinewimage --size "$SIZE" -f 4.0 T1_csf.nii
qinewimage --size "$SIZE" -f 2.5 T2_csf.nii
qinewimage --size "$SIZE" -f "0.18" tau_m.nii
qinewimage --size "$SIZE" -g "0 0. 200." f0.nii
qinewimage --size "$SIZE" -g "1 0.75 1.25" B1.nii
qinewimage --size "$SIZE" -g "2 0.05 0.25" f_m.nii
qinewimage --size "$SIZE" -f 0.05 f_csf.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="3 4 5 6 7 9 13 18"
SPGR_TR="0.0065"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 16 21 27 33 40 51 68"
SSFP_PC="180 0"
SSFP_TR="0.005"

INPUT="$SPGR_FILE
SPGR
$SPGR_FLIP
$SPGR_TR
$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END
"

echo "$INPUT" > mcd_input.txt
NOISE="0.002"
qisignal --model=3 -v --noise=$NOISE<< END_SIG
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
END_SIG

qimcdespot $OPTS -M3 -bB1.nii -ff0.nii -r -v -n < mcd_input.txt
qidiff --baseline=f_m.nii --input=3C_f_m.nii --noise=$NOISE --tolerance=250 --verbose

}