# Copyright Tobias Wood 2017
# Tests for phase unwrapping

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Path-based Phase Unwrapping" {

SIZE="32,32,32"
# Phase unwrapping is not unique circa global multiples of 2*pi
# The following ramp & wrapped values produce maps the coincide with the current version
qinewimage --size=$SIZE --grad="0 -6.2832 18.8496" ramp$EXT
qinewimage --size=$SIZE --grad="0 -12.566 12.566" --wrap=6.283 wrapped$EXT
qi_unwrap_path wrapped$EXT --verbose
qidiff --baseline=ramp$EXT --input=wrapped_unwrapped$EXT --tolerance=1 --abs --verbose

}

# @test "Laplacian Phase Unwrapping" {

# SIZE="32,32,32"
# qinewimage --size=$SIZE --grad="0 0 12.566" ramp$EXT
# qinewimage --size=$SIZE --grad="0 0 12.566" --wrap=6.283 wrapped$EXT
# qi_unwrap_laplace wrapped$EXT --verbose
# qidiff --baseline=ramp$EXT --input=wrapped_unwrapped$EXT --tolerance=1 --abs --verbose

# }