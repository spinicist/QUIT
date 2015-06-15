# QUIT - QUantitative Imaging Tools #

Credit / Blame / Contact - Tobias Wood - tobias.wood@kcl.ac.uk

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
If you find the tools useful the author would love to hear from you.

This is the updated version of QUIT based on ITK http://www.itk.org
The previous version is here http://github.com/spinicist/old_quit

# Brief Description #

A collection of programs for processing quantitative MRI data, in particular
DESPOT and mcDESPOT.

# Thanks #

Many thanks to Samuel Hurley for many, many suggestions.
Thanks to Anna Combes, Amy McDowell and Sjoerd Vos for finding bugs.

# Installation #

## Dependencies & Requirements ##

1. A C++11 compliant compiler (GCC 4.8.0+, Clang 3.1+)
2. CMake version 3.2 or higher (http://www.cmake.org)
3. Eigen version 3.2.4 or higher (http://eigen.tuxfamily.org)
4. ITK version 4.7.2 (http://www.itk.org)

WARNING - You will require a recent compiler for DESPOT as it uses C++11
features. GCC 4.8.0 is a MINIMUM as earlier versions do not support C++11
threads. Clang 3.1 will also work. You can find out the version of your
system gcc by running "gcc -v".

## Compilation ##

There are two suggested methods to build the tools. If you are unfamiliar with
CMake and command-line building use method 1. If you are familiar with them
then you can use method 2.

1. Download and untar/unzip QUIT. In a terminal window, change directory to
where you unpacked QUIT. Then type './easy_build.sh'. This will download Eigen
& ITK into subdirectories, and then compile QUIT in a '/Build' subdirectory.
2. Download and compile ITK & Eigen. Download and unpack QUIT. Create a build
directory, e.g. 'QUIT/Build'. Then run CMake as standard, e.g. 'ccmake ../'
from within 'QUIT/Build'. Specify the directories of ITK and Eigen. It is
likely that CMake will report an error on the first 'Configure' attempt if it
cannot locate Eigen. After configuring and generating, run 'make -j'.

# Usage #

Each product has some basic usage instructions that will be printed with either
the -h or --help options, e.g. “qimcdespot -h". The programs will then prompt you
to enter flip-angles, TRs, etc. as they require based on your input options.

Once you are used to the input format, I recommend writing down the input in a
plain text file, and then using the redirect operator (“<“) to pass the input
to your program. If you do this, you can suppress the prompts to enter input
with the “-n” option. So to run DESPOT1-HIFI you would type:

	qidespot1hifi -m mask_file.nii -n spgr_file.nii irspgr_file.nii < in.txt

where “in.txt” is a file containing:

	# A # symbol means a comment
	3 18   # SPGR flip-angles
	0.0083 # SPGR TR in SECONDS, not milli-seconds
	5      # IR-SPGR readout flip-angle
	0.0083 # IR-SPGR readout TR
	92     # IR-SPGR echo train length
	0.45   # IR-SPGR TI in SECONDS



