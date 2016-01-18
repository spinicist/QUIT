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
int print_all = true;
int print_size = false, print_spacing = false, print_origin = false, print_direction = false;

const string usage {
"Usage is: qihdr [options] input[s] \n\
\n\
By default, a summary of the header is printed. If options below are specified,\n\
only those parts of the header will be printed. Multiple files can be input,\n\
in which case the header info is written for each in order.\n\
\n\
Options:\n\
    --help, -h  : Print this message.\n\
    --size      : Print the image dimensions.\n\
    --spacing   : Print the image spacing.\n\
    --origin    : Print the origin.\n\
    --direction : Print the image direction.\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"size", no_argument, 0, 's'},
    {"spacing", no_argument, 0, 'p'},
    {"origin", no_argument, 0, 'o'},
    {"direction", no_argument, 0, 'd'},
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
            case 's': print_size = true; print_all = false; break;
            case 'p': print_spacing = true; print_all = false; break;
            case 'o': print_origin = true; print_all = false; break;
            case 'd': print_direction = true; print_all = false; break;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }

    while (optind < argc) {
        QI::ReadImageF::Pointer thisImg = QI::ReadImageF::New();
        if (print_all) cout << "Header for: " << string(argv[optind]) << endl;
        thisImg->SetFileName(argv[optind++]);
        thisImg->Update();
        if (print_all) cout << "Size:      "; if (print_all || print_size)      cout << thisImg->GetOutput()->GetLargestPossibleRegion().GetSize() << endl;
        if (print_all) cout << "Spacing:   "; if (print_all || print_spacing)   cout << thisImg->GetOutput()->GetSpacing() << endl;
        if (print_all) cout << "Origin:    "; if (print_all || print_origin)    cout << thisImg->GetOutput()->GetOrigin() << endl;
        if (print_all) cout << "Direction: " << endl; if (print_all || print_direction) cout << thisImg->GetOutput()->GetDirection() << endl;
    }
    return EXIT_SUCCESS;
}
