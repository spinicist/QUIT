# Copyright Tobias Wood 2017
# Tests for perfusion (ASL)

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Perfusion (ASL)" {

SIZE="32, 32, 32, 2"
qinewimage --verbose --dims=4 --size="$SIZE" --step="3 0 1 2" asl.nii
qi_asl --verbose asl.nii <<END_INPUT
{
    "sequence" : {
        "TR": 4.0,
        "label_time": 3.0,
        "post_label_delay": 0.3
    }
}
END_INPUT

}
