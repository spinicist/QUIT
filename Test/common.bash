#DATADIR="$BATS_TEST_DIRNAME/basic"
DATADIR="$BATS_TMPDIR"
BINDIR="$BATS_TEST_DIRNAME/.."

init_tests() {
    mkdir -p $DATADIR
    cd $DATADIR
    echo "Running tests in $DATADIR"
    echo "Removing any previous image files"
    if ls *.nii* > /dev/null 2>&1; then
        rm *.nii*
    fi
}

img_diff() {
    INPUT="$1"
    REFERENCE="$2"
    ADDED_NOISE="$3"
    TOLERANCE="$4"
    fslmaths $INPUT -sub $REFERENCE -div $REFERENCE -sqr diff_sqr
    MEAN_SQR=$( fslstats diff_sqr -m )
    REL_NOISE=$( echo $MEAN_SQR $ADDED_NOISE | awk ' { print sqrt($1)/$2 } ' )
    # Check for nan/inf/etc. because on some platforms awk will treat these as 0
    REGEXP='^-?[0-9]+([.][0-9]+)?(e-?[0-9]+)?$'
    if ! [[ $REL_NOISE =~ $REGEXP ]] ; then
        echo "Comparison test $NAME failed, mean diff $REL_NOISE is not a valid number"
        return 1
    fi
    # Now do the tolerance test
    TEST=$(echo "$REL_NOISE $TOLERANCE" | awk ' { if(($1)<=($2)) { print 1 } else { print 0 }}')
    if [ "$TEST" -eq "1" ]; then
        echo "Images match within tolerance ($REL_NOISE < $TOLERANCE)"
        return 0
    else
        echo "Images do not match within tolerance ($REL_NOISE > $TOLERANCE)"
        return 1
    fi
}