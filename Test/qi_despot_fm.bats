# Copyright Tobias Wood 2017
# Tests for DESPOT2-FM

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "DESPOT2-FM" {

SIZE="16,16,9"
qinewimage PD$EXT --size "$SIZE" -f 1.0
qinewimage T1$EXT --size "$SIZE" -f 1.0
qinewimage T2$EXT --size "$SIZE" -g "1 0.025 0.125"
qinewimage f0$EXT --size "$SIZE" -g "2 -150.0 150.0"
qinewimage B1$EXT --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp$EXT"
SSFP_FLIP="12,12,12,12,65,65,65,65"
SSFP_PINC="0,90,180,270,0,90,180,270"
SSFP_TR="0.005"
SSFP_Trf="0.001"

NOISE="0.002"
qisignal --verbose --model=1 --noise=$NOISE $SSFP_FILE << END_SIG
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "T2$EXT",
    "f0": "f0$EXT",
    "B1": "",
    "SequenceGroup": {
        "sequences": [
            {
                "SSFP": {
                    "TR": $SSFP_TR,
                    "FA": [$SSFP_FLIP],
                    "PhaseInc": [$SSFP_PINC]
                }
            }
        ]
    }
}
END_SIG

qidespot2fm --verbose -bB1.nii T1.nii ${SSFP_FILE} --asym << END_FM
{
    "SSFP": {
        "TR": $SSFP_TR,
        "PhaseInc": [$SSFP_PINC],
        "FA": [$SSFP_FLIP]
    }
}
END_FM

qidiff --baseline=T2.nii --input=FM_T2.nii --noise=$NOISE --tolerance=50 --verbose

}
