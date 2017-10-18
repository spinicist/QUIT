# Copyright Tobias Wood 2017
# Tests for polynomial fitting

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Polynomial Fitting" {

SIZE="32, 32, 32"
qinewimage --size "$SIZE" --step="2 0 1 2" mask.nii
qipolyimg mask.nii baseline.nii --order=2 -v << END_INPUT
0 0 0
0 1 2 4 1 0 0 1 0 1
END_INPUT
qipolyfit baseline.nii --order=2 --print-terms > unmasked_terms.txt
qipolyfit baseline.nii --order=2 --print-terms --mask=mask.nii > masked_terms.txt
qipolyfit baseline.nii --order=2 --print-terms --robust > robust_terms.txt
qipolyimg mask.nii unmasked.nii --order=2 < unmasked_terms.txt
qipolyimg mask.nii masked.nii --order=2 < masked_terms.txt
qipolyimg mask.nii robust.nii --order=2 < robust_terms.txt
img_diff baseline.nii unmasked.nii 1 1
img_diff baseline.nii masked.nii   1 1
img_diff baseline.nii unmasked.nii 1 1

}
