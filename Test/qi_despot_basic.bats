# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "DESPOT1-Basic" {

# Setup parameters
SPGR_FILE="spgr$EXT"
SPGR_FLIP="3,3,20,20"
SPGR_TR="0.01"
SIZE="16,16,16"
NOISE="0.01"
qinewimage --size "$SIZE" -g "1 0.8 1.0" PD$EXT
qinewimage --size "$SIZE" -g "0 0.5 1.5" T1$EXT
qisignal --model=1 -v --noise=$NOISE $SPGR_FILE << OUT
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "",
    "f0": "",
    "B1": "",
    "Sequences": [
        {
            "SPGR": {
                "TR": $SPGR_TR,
                "FA": [$SPGR_FLIP]
            }
        }
    ]
}
OUT
qidespot1 $SPGR_FILE --verbose <<OUT
{
    "SPGR": {
        "TR": $SPGR_TR,
        "FA": [$SPGR_FLIP]
    }
}
OUT
qidiff --baseline=T1.nii --input=D1_T1.nii --noise=$NOISE --tolerance=35 --verbose

}

@test "DESPOT2-Basic" {

# Setup parameters
SSFP_FILE="ssfp$EXT"
SSFP_FLIP="15,60"
SSFP_PINC="180,180"
SSFP_TR="0.01"
SIZE="16,16,16"
NOISE="0.01"
qinewimage --size "$SIZE" -g "1 0.8 1.0" PD.nii
[ -e PD.nii ]
qinewimage --size "$SIZE" -g "0 0.5 1.5" T1.nii
[ -e T1.nii ]
qinewimage --size "$SIZE" -g "2 0.02 0.1" T2.nii
[ -e T2.nii ]
qisignal --model=1 -v --noise=$NOISE $SSFP_FILE << OUT
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "T2$EXT",
    "f0": "",
    "B1": "",
    "Sequences": [
        {
            "SSFP": {
                "TR": $SSFP_TR,
                "FA": [$SSFP_FLIP],
                "PhaseInc": [$SSFP_PINC]
            }
        }
    ]
}
OUT
qidespot2 T1.nii $SSFP_FILE --verbose <<OUT
{
    "SSFP": {
        "TR": $SSFP_TR,
        "PhaseInc": [$SSFP_PINC],
        "FA": [$SSFP_FLIP]
    }
}
OUT
qidiff --baseline=T2.nii --input=D2_T2.nii --noise=$NOISE --tolerance=30 --verbose

}