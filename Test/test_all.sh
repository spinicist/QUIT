#!/bin/bash

# Tobias Wood 2015
# Driver for all tests in order

source ./test_common.sh

run_test "SCDESPOT" ./test_scdespot.sh
run_test "MCDESPOT" ./test_mcdespot.sh
run_test "RELAX" ./test_relax.sh

echo "Finished tests."
