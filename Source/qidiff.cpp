/*
 *  qidiff.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "QI/Util.h"
#include "QI/Args.h"

using namespace std;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qidiff file1 file2 [options]\n\
\n\
A tool to compare images, intended for use with the library tests.\n\
Subtracts one image from the other, checks for NaNs etc., then\n\
calculates the root sum-of-squares residuals and compares to a\n\
threshold (default 0)."
};
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {

    return EXIT_SUCCESS;
}


