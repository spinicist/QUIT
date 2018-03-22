# Copyright Tobias Wood 2018
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "SSFP Mag Banding" {
    SIZE="64,64,16"
    qinewimage --size "$SIZE" -f 1.0 PD$EXT
    qinewimage --size "$SIZE" -g "1 0.5 1.5" T1$EXT
    qinewimage --size "$SIZE" -g "2 0.04 0.1" T2$EXT
    qinewimage --size "$SIZE" -g "0 -200 200" f0$EXT
    qisignal --model=1 -v --noise=0.01 ssfp_mag$EXT << OUT
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
                    "TR": 0.01,
                    "FA": [30,30,30,30],
                    "PhaseInc": [180,0,90,270]
                }
            }
        ]
    }
}
OUT
    qi_ssfp_bands ssfp_mag$EXT --method=M --magnitude --verbose
    qi_ssfp_bands ssfp_mag$EXT --method=N --magnitude --verbose
    # No proper ground truth here, compare the methods to each other
    qidiff --baseline=ssfp_mag_MagMean$EXT --input=ssfp_mag_Max$EXT --noise=0.01 --tolerance=30 --verbose
}

@test "SSFP GS Banding" {
    SIZE="64,64,16"
    qinewimage --size "$SIZE" -f 1.0 PD$EXT
    qinewimage --size "$SIZE" -g "1 0.5 1.5" T1$EXT
    qinewimage --size "$SIZE" -g "2 0.04 0.1" T2$EXT
    qinewimage --size "$SIZE" -g "0 -200 200" f0$EXT
    qisignal --model=1 -v --noise=0.01 ssfp$EXT ssfpgs$EXT --complex << OUT
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
                    "TR": 0.01,
                    "FA": [30,30,30,30],
                    "PhaseInc": [180,0,90,270]
                }
            },{
                "SSFPGS": {
                    "TR": 0.01,
                    "FA": [30]
                }
            }
        ]
    }
}
OUT

    qi_ssfp_bands ssfp$EXT --method=G --magnitude --alt-order --verbose
    qi_ssfp_bands ssfp$EXT --method=G -2 --magnitude --alt-order --verbose
    qicomplex -x ssfpgs$EXT -M ssfpgs_mag$EXT
    qidiff --baseline=ssfpgs_mag$EXT --input=ssfp_GSL$EXT --noise=0.01 --tolerance=15 --verbose
    qidiff --baseline=ssfpgs_mag$EXT --input=ssfp_GSL2$EXT --noise=0.01 --tolerance=13 --verbose
}