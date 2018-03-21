# Copyright Tobias Wood 2018
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "SSFP Planet" {
    SIZE="64,64,16"
    qinewimage --size "$SIZE" -f 1.0 PD$EXT
    qinewimage --size "$SIZE" -g "1 0.5 1.5" T1$EXT
    qinewimage --size "$SIZE" -g "2 0.05 0.1" T2$EXT
    qinewimage --size "$SIZE" -g "0 -200 200" f0$EXT
    qisignal --model=1 -v --noise=0.00 ssfp$EXT --complex << INPUT
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "T2$EXT",
    "f0": "",
    "B1": "",
    "SequenceGroup": {
        "sequences": [
            {
                "SSFPEcho": {
                    "TR": 0.01,
                    "FA": [30,30,30,30,30,30],
                    "PhaseInc": [180,240,300,0,60,120]
                }
            }
        ]
    }
}
INPUT
    qi_ssfp_ellipse ssfp$EXT --verbose <<INPUT
{
    "SSFPEllipse": {
        "TR": 0.01,
        "FA": [30],
        "PhaseInc": [180,240,300,0,60,120]
    }
}
INPUT
    qi_ssfp_planet ES_G$EXT ES_a$EXT ES_b$EXT --verbose << INPUT
{
    "SSFPGS": {
        "TR": 0.01,
        "FA": [30]
    }
}
INPUT
    # No proper ground truth here, compare the methods to each other
    qidiff --baseline=T1$EXT --input=PLANET_T1$EXT --noise=0.01 --tolerance=30 --verbose
    qidiff --baseline=T2$EXT --input=PLANET_T2$EXT --noise=0.01 --tolerance=30 --verbose
    qidiff --baseline=PD$EXT --input=PLANET_PD$EXT --noise=0.01 --tolerance=30 --verbose
}
