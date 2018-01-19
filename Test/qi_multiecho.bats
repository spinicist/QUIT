# Copyright Tobias Wood 2017
# Tests for relaxometry programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Multi-Echo" {

SIZE="16,16,31"
qinewimage PD.nii --size "$SIZE" -f 1
qinewimage T1.nii --size "$SIZE" -g "0 0.5 1.5"
qinewimage T2.nii --size "$SIZE" -g "1 0.01 0.1"
qinewimage f0.nii --size "$SIZE" -g "2 -25.0 25.0"
qinewimage B1.nii --size "$SIZE" -f 1

# Setup parameters
SPIN_FILE="me.nii"
SPIN_PAR="2.5
0.005
0.005
16"

# Create input for Single Component
echo "$MCSIG_INPUT" > qisignal.in
NOISE="0.002"
qisignal --model=1 -n -v --noise=$NOISE << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SPIN_FILE
SPINECHO
$SPIN_PAR
END
END_SIG

echo "$SPIN_PAR" > multiecho.in
qimultiecho $SPIN_FILE -n -v -al -oLL_ < multiecho.in
qimultiecho $SPIN_FILE -n -v -an -oLM_ < multiecho.in
qimultiecho $SPIN_FILE -n -v -aa -oAR_   < multiecho.in

qidiff --baseline=T2.nii --input=LL_ME_T2.nii --noise=$NOISE --tolerance=50 --verbose
qidiff --baseline=T2.nii --input=LM_ME_T2.nii --noise=$NOISE --tolerance=50 --verbose
qidiff --baseline=T2.nii --input=AR_ME_T2.nii --noise=$NOISE --tolerance=50 --verbose

}
