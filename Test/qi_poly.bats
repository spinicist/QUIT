# Copyright Tobias Wood 2017
# Tests for polynomial fitting

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Polynomial Fitting" {

SIZE="32, 32, 32"
qinewimage --size "$SIZE" --step="2 0 1 2" mask.nii
# Set constant term to 1 so there are no zero voxels, because qidiff doesn't like that
qipolyimg mask.nii baseline.nii --order=2 -v << END_INPUT
{
    "center": [0, 0, 0],
    "scale": 1,
    "coeffs": [1, 1, 2, 4, 1, 0, 0, 1, 0, 1]
}
END_INPUT
qipolyfit baseline.nii --order=2 --print-terms > unmasked_terms.txt
qipolyfit baseline.nii --order=2 --print-terms --mask=mask.nii > masked_terms.txt
qipolyfit baseline.nii --order=2 --print-terms --robust > robust_terms.txt
qipolyimg mask.nii unmasked.nii --order=2 < unmasked_terms.txt
qipolyimg mask.nii masked.nii --order=2 < masked_terms.txt
qipolyimg mask.nii robust.nii --order=2 < robust_terms.txt
qidiff --baseline=baseline.nii --input=unmasked.nii --tolerance=1 --verbose
qidiff --baseline=baseline.nii --input=masked.nii   --tolerance=1 --verbose
qidiff --baseline=baseline.nii --input=unmasked.nii --tolerance=1 --verbose

}
