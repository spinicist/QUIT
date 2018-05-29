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

SPIN_SEQUENCE="
\"MultiEcho\" : {
    \"TR\" : 2.5,
    \"TE1\" : 0.005,
    \"ESP\" : 0.005,
    \"ETL\" : 16
}
"

SPIN_FILE="me.nii"
NOISE="0.02"
qisignal --model=1 -v --noise=$NOISE $SPIN_FILE << END_SIG
{
    "PD" : "PD.nii",
    "T1" : "T1.nii",
    "T2" : "T2.nii",
    "f0" : "f0.nii",
    "B1" : "B1.nii",
    "SequenceGroup" : {
        "sequences" : [
            {
                $SPIN_SEQUENCE
            }
        ]
    }
}
END_SIG

echo { "$SPIN_SEQUENCE" } > multiecho.in
qimultiecho $SPIN_FILE -v -al -oLL_ < multiecho.in
qimultiecho $SPIN_FILE -v -an -oLM_ < multiecho.in
qimultiecho $SPIN_FILE -v -aa -oAR_   < multiecho.in

qidiff --baseline=T2.nii --input=LL_ME_T2.nii --noise=$NOISE --tolerance=20 --verbose
qidiff --baseline=T2.nii --input=LM_ME_T2.nii --noise=$NOISE --tolerance=5 --verbose
qidiff --baseline=T2.nii --input=AR_ME_T2.nii --noise=$NOISE --tolerance=5 --verbose

}
