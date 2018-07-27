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
