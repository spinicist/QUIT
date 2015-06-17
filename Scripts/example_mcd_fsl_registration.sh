#!/bin/bash -e

#
# Example Processing Script For MCDESPOT data
# By Tobias Wood, with help from Anna Coombes and Sam Hurley
#

if [ $# -neq 0 ]; then
echo << END_USAGE
Usage: $0 file1.nii file2.nii [file3.nii ...]

This script will register all volumes in each file to the first volume in the
first file specified.

Requires FSL.

This is almost certainly not the optimal way to do the registration, but should
be regarded as a framework to improve on.

END_USAGE
exit 1;
fi

TARGET="target"

#
# Extract a target volume
#

fslroi $1 $TARGET 0 1

while test $# -gt 0
do
	FILE="$1"
	echo "Splitting $FILE into volumes..."

	fslsplit $FILE ${FILE%%.nii*}_split -t

	echo "Registering volumes"
	for VOL in ${FILE%%.nii*}_split*; do
		flirt -dof 6 -cost mutualinfo -in $VOL -ref $TARGET -out ${VOL%%nii*}_reg
	done

	echo "Merging into registered file"
	fslmerge -t ${FILE%%.nii*}_reg ${VOL$$.nii*}_reg*

	echo "Removing temporary files"
	rm ${FILE%%.nii*}_split*
	shift # Prep for next iteration
done
rm $TARGET
