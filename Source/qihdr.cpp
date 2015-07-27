/*
 *  qihdr.cpp
 *
 *  Created by Tobias Wood on 11/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>

#include "Types.h"
#include "Util.h"

using namespace std;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qihdr [options] input[s] \n\
\n\
Options:\n\
    --help, -h        : Print this message.\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};
const char *short_options = "h";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'h':
                cout << usage << endl;
                return EXIT_SUCCESS;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }

    while (optind < argc) {
        QI::ReadImageF::Pointer thisImg = QI::ReadImageF::New();
        cout << "Header for: " << string(argv[optind]) << endl;
        thisImg->SetFileName(argv[optind++]);
        thisImg->Update();
        cout << "Size:      " << thisImg->GetOutput()->GetLargestPossibleRegion().GetSize() << endl;
        cout << "Spacing:   " << thisImg->GetOutput()->GetSpacing() << endl;
        cout << "Origin:    " << thisImg->GetOutput()->GetOrigin() << endl;
        cout << "Direction: " << thisImg->GetOutput()->GetDirection() << endl;
    }
    return EXIT_SUCCESS;
}
