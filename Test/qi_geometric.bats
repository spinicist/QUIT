# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Geometric Solution" {

# First, create input data
SIZE="16,16,51"
qinewimage PD$EXT --size "$SIZE" -f "1.0"
qinewimage T1$EXT --size "$SIZE" -g "0 1.0 2.0"
qinewimage T2$EXT --size "$SIZE" -g "1 0.025 0.25"
qinewimage f0$EXT --size "$SIZE" -g "2 -200.0 200.0"
qinewimage B1$EXT --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp$EXT"
GS_FILE="gs_ref$EXT"
SSFP_TR="0.005"
SSFP_FLIP="45,45,45,45"
SSFP_PINC="180,270,0,90"
NOISE="2e-3"
qisignal --verbose -x --noise=$NOISE -T1 $SSFP_FILE $GS_FILE << END_SIG
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
            },
            {
                "SSFPGS": {
                    "TR": $SSFP_TR,
                    "FA": [$SSFP_FLIP],
                    "PhaseInc": [$SSFP_PINC]
                }
            }
        ]
    }
}
END_SIG

qissfpbands ${SSFP_FILE} --ph-incs=4 --verbose
[ -e ssfp_GSL$EXT ]
qicomplex -x ssfp_GSL$EXT -M gs_mag$EXT -P gs_ph$EXT --verbose
[ -e gs_mag$EXT ]
[ -e gs_ph$EXT ]
qicomplex -x gs_ref$EXT -M gs_ref_mag$EXT --verbose
[ -e gs_ref_mag$EXT ]
qidiff --baseline=gs_ref_mag$EXT --input=gs_mag$EXT --noise=$NOISE --tolerance=50 --verbose

}
