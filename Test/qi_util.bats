# Copyright Tobias Wood 2018
# Tests for various utilities

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "RF Profile" {

SIZE="32, 32, 32"
qinewimage --size "$SIZE" --grad="1 0 1" grad.nii
qi_rfprofile grad.nii b1.nii <<INPUT
{
    "rf_pos" :  [0, 1],
    "rf_vals" : [ [0, 1],
                  [1, 0] ]
}
INPUT
}
