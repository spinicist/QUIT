# Copyright Tobias Wood 2017
# Tests for polynomial fitting

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Polynomial Fitting" {

SIZE="32, 32, 32"
qinewimage --size "$SIZE" --grad_dim=2 --grad_vals=0,1 --steps=2 mask.nii
# Set constant term to 1 so there are no zero voxels, because qidiff doesn't like that
qipolyimg mask.nii baseline.nii --order=2 -v << END_INPUT
{
    "center": [0, 0, 0],
    "scale": 1,
    "coeffs": [1, 1, 2, 4, 1, 0, 0, 1, 0, 1]
}
END_INPUT
qipolyfit baseline.nii --order=2 --print-terms > unmasked_terms.json
qipolyfit baseline.nii --order=2 --print-terms --mask=mask.nii > masked_terms.json
qipolyfit baseline.nii --order=2 --print-terms --robust > robust_terms.json
qipolyimg mask.nii unmasked.nii --order=2 < unmasked_terms.json
qipolyimg mask.nii masked.nii --order=2 < masked_terms.json
qipolyimg mask.nii robust.nii --order=2 < robust_terms.json
qidiff --baseline=baseline.nii --input=unmasked.nii --verbose
qidiff --baseline=baseline.nii --input=masked.nii   --verbose
qidiff --baseline=baseline.nii --input=unmasked.nii --verbose

}
