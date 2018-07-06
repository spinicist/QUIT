# Copyright Tobias Wood 2017
# Tests for perfusion (ASL)

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Perfusion (ASL)" {

SIZE="32,32,32,2"
qinewimage --verbose --dims=4 --size="$SIZE" --step="3 1 1.06 2" asl.nii
qi_asl --verbose asl.nii <<END_INPUT
{
    "CASL" : {
        "TR": 4.0,
        "label_time": 3.0,
        "post_label_delay": [0.3]
    }
}
END_INPUT
}

@test "Perfusion (ASE)" {
# Currently, this test uses a T2 decay curve as the input and checks that R2' is equal to 1/T2
SIZE="16,16,16"
qinewimage PD.nii --size "$SIZE" -f 1
qinewimage T1.nii --size "$SIZE" -f 1
qinewimage T2.nii --size "$SIZE" -f 0.1
qinewimage R2.nii --size "$SIZE" -f 10

SPIN_SEQUENCE="
\"MultiEchoFlex\" : {
    \"TR\" : 2.5,
    \"TE\" : [ 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07 ]
}
"
SPIN_FILE="me.nii"
NOISE="0.02"
qisignal --model=1 -v --noise=$NOISE $SPIN_FILE << END_SIG
{
    "PD" : "PD.nii",
    "T1" : "T1.nii",
    "T2" : "T2.nii",
    "f0" : "",
    "B1" : "",
    "Sequences" : [
        { $SPIN_SEQUENCE }
    ]
}
END_SIG

qi_ase_oef --verbose $SPIN_FILE << END_INPUT
{
    $SPIN_SEQUENCE
}
END_INPUT
qidiff --baseline=R2.nii --input=me_R2prime.nii --noise=$NOISE --tolerance=5 --verbose

}
