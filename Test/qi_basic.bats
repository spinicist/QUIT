# Copyright Tobias Wood 2017
# Test scripts for the basic QUIT Tools

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Create Image" {
    # Setup parameters
    command echo $DATADIR
    SIZE="32,32,32"
    ONE_SIZE="32"
    IMAGE="test_image$EXT"
    qinewimage --size "$SIZE" "$IMAGE"
    [ -e "$IMAGE" ]
    RESULT="$( qihdr "$IMAGE" --size=1 )"
    [ "$RESULT"  -eq "$ONE_SIZE" ]
}

@test "Simulate Sequence" {
    # Setup parameters
    command echo $DATADIR
    SIZE="8,8,8"
    qinewimage --size "$SIZE" PD$EXT
    qinewimage --size "$SIZE" T1$EXT
    qisignal --model=1 -v --noise=$NOISE SPGR$EXT << OUT
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "",
    "f0": "",
    "B1": "",
    "Sequences": [
        {
            "SPGR": {
                "TR": 0.01,
                "FA": [5, 10]
            }
        }
    ]
}
OUT
    [ -e "SPGR$EXT" ]
}
