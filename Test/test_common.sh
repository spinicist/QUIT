#!/bin/bash -e

# Tobias Wood 2015
# Common functions for QUITK Tests

# Check for presence of FSL
if [ -n "$(which fslmaths)" ]; then
	HAVE_FSL="1"
	echo "FSL detected."
	export FSLOUTPUTTYPE=NIFTI
else
	HAVE_FSL="0"
	echo "FSL not detected. Cannot run comparison tests."
fi

# Simple test function
SILENCE_TESTS=0
function run_test {
	# $1 is test name, remainder is command to run
	NAME="$1"
	shift
	printf "Starting test $NAME..."
	if [ "$SILENCE_TESTS" -eq "1" ]; then
		"$@" > "$NAME.log"
	else
		"$@"
	fi
	local STATUS=$?
	if [ $STATUS -ne 0 ]; then
		printf "Failed.\n" >&2
		exit $STATUS
	else
		printf "Passed.\n" >&1
	fi
	return $STATUS
}

# Use FSL to compare 2 images and check if the average difference is below the tolerance
function compare_test {
	# $1 test name
	# $2 reference, $3 test image, $4 tolerance
	NAME="$1"
	REF=$2
	TEST=$3
	TOL=$4
	DIFF=${REF%.nii}_${TEST%.nii}
	if [ "$HAVE_FSL" -eq "1" ]; then
		fslmaths $REF -sub $TEST $DIFF
		MEAN=$( fslstats $DIFF -M )
		STD=$( fslstats $DIFF -s )
		ABSMEAN=$(echo $MEAN | awk ' { print sqrt($1*$1) } ' )
		# Check for nan/inf/etc. because on some platforms awk will treat these as 0
		REGEXP='^-?[0-9]+([.][0-9]+)?(e-?[0-9]+)?$'
		if ! [[ $ABSMEAN =~ $REGEXP ]] ; then
			echo "Comparison test $NAME failed, mean diff $ABSMEAN is not a valid number"
			exit 1
		fi
		# Now do the tolerance test
		TEST=$(echo "$ABSMEAN $TOL" | awk ' { if(($1)<=($2)) { print 1 } else { print 0 }}')
		if [ "$TEST" -eq "1" ]; then
			echo "Comparison test $NAME passed, mean diff was $ABSMEAN tolerance $TOL (std was $STD)"
		else
			echo "Comparison test $NAME failed, mean diff was $ABSMEAN tolerance $TOL (std was $STD)"
			exit 1
		fi
	else
		echo "FSL not present, skipping test $NAME"
	fi
	return 0
}

# Setup environment
QUITDIR=$PWD/..
export QUIT_EXT=NIFTI
#QUITVER=$(cat ../src/version | sed -e 's/^"//'  -e 's/"$//')
