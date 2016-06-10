#!/bin/bash -eu

# Tobias Wood 2015
# Driver for all tests in order

source ./test_common.sh

run_test "DESPOT1" ./test_d1.sh
run_test "FM"      ./test_fm.sh
run_test "SSFP"    ./test_ssfp.sh
run_test "RELAX"   ./test_relax.sh
run_test "2C"      ./test_2C.sh
run_test "3C"      ./test_3C.sh

echo "Finished tests."
