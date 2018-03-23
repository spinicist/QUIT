# Copyright Tobias Wood 2018
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "SSFP PLANET" {
    SIZE="64,64,16"
    NOISE="0.005"
    qinewimage --size "$SIZE" -f 1.0 PD$EXT
    qinewimage --size "$SIZE" -g "0 0.5 1.5" T1$EXT
    qinewimage --size "$SIZE" -g "1 0.05 0.1" T2$EXT
    qinewimage --size "$SIZE" -g "2 -200 200" f0$EXT
    qisignal --model=1 -v --noise=$NOISE ssfp$EXT --complex << INPUT
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "T2$EXT",
    "f0": "f0$EXT",
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
    qi_ssfp_ellipse ssfp$EXT --verbose --algo=d << INPUT
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
    qidiff --baseline=T1$EXT --input=PLANET_T1$EXT --noise=$NOISE --tolerance=40 --verbose
    qidiff --baseline=T2$EXT --input=PLANET_T2$EXT --noise=$NOISE --tolerance=35 --verbose
    qidiff --baseline=PD$EXT --input=PLANET_PD$EXT --noise=$NOISE --tolerance=20 --verbose
}

@test "SSFP EMT" {
    SIZE="64,64,16"
    NOISE="0.005"
    qinewimage --size "$SIZE" -f 1.0 PD$EXT
    qinewimage --size "$SIZE" -g "0 0.5 1.5" T1$EXT
    qinewimage --size "$SIZE" -g "1 0.05 0.1" T2$EXT
    qinewimage --size "$SIZE" -g "2 -200 200" f0$EXT
    qisignal --model=1 -v --noise=$NOISE ssfp$EXT --complex << INPUT
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "T2$EXT",
    "f0": "f0$EXT",
    "B1": "",
    "SequenceGroup": {
        "sequences": [
            {
                "SSFPEcho": {
                    "TR": 0.01,
                    "FA": [16,16,16,16,16,16,32,32,32,32,32,32],
                    "PhaseInc": [180,240,300,0,60,120,180,240,300,0,60,120]
                }
            }
        ]
    }
}
INPUT
    qi_ssfp_ellipse ssfp$EXT --verbose --algo=d << INPUT
{
    "SSFPEllipse": {
        "TR": 0.01,
        "FA": [16,32],
        "PhaseInc": [180,240,300,0,60,120]
    }
}
INPUT
    qi_ssfp_emt ES_G$EXT ES_a$EXT ES_b$EXT --verbose << INPUT
{
    "SSFPMT": {
        "TR": [0.01, 0.01],
        "Trf": [0.001, 0.001],
        "FA": [16, 32],
        "intB1": [0.005904, 0.005904],
        "PhaseInc": [180, 240, 300, 0, 60, 120]
    }
}
INPUT
    # No proper ground truth here, compare the methods to each other
    qidiff --baseline=T1$EXT --input=EMT_T1f$EXT --noise=$NOISE --tolerance=40 --verbose
    qidiff --baseline=T2$EXT --input=EMT_T2f$EXT --noise=$NOISE --tolerance=35 --verbose
    qidiff --baseline=PD$EXT --input=EMT_M0$EXT --noise=$NOISE --tolerance=20 --verbose
}
