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

#include <getopt.h>
#include <iostream>

#include "QI/Types.h"
#include "QI/Util.h"

using namespace std;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
int print_all = true;
int print_dims = false, print_type = false, print_size = false, print_spacing = false, print_origin = false, print_direction = false;
int dim3 = false;
const string usage {
"Usage is: qihdr [options] input[s] \n\
\n\
By default, a summary of the header is printed. If options below are specified,\n\
only those parts of the header will be printed. Multiple files can be input,\n\
in which case the header info is written for each in order.\n\
\n\
Options:\n\
    --help, -h  : Print this message.\n\
    --3D        : Discard info for higher dimensions.\n\
    --dims      : Print the number of dimensions.\n\
    --dtype     : Print the data type.\n\
    --size      : Print the image dimension sizes.\n\
    --spacing   : Print the image spacing.\n\
    --origin    : Print the origin.\n\
    --direction : Print the image direction.\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"3D", no_argument, &dim3, true},
    {"dims", no_argument, 0, 'D'},
    {"dtype", no_argument, 0, 'T'},
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
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case 'D': print_dims = true; print_all = false; break;
            case 'T': print_type = true; print_all = false; break;
            case 's': print_size = true; print_all = false; break;
            case 'p': print_spacing = true; print_all = false; break;
            case 'o': print_origin = true; print_all = false; break;
            case 'd': print_direction = true; print_all = false; break;
            case 0: break; // A flag was set
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }

    while (optind < argc) {
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(argv[optind], itk::ImageIOFactory::ReadMode);
        if (!imageIO) {
            cerr << "Could not open: " << string(argv[optind]) << endl;
            break;
        }
        imageIO->SetFileName(string(argv[optind]));
        imageIO->ReadImageInformation();
        size_t dims = imageIO->GetNumberOfDimensions();
        if (print_all) cout << "File:        " << string(argv[optind]) << endl;
        if (print_all) cout << "Dimension:   "; if (print_all | print_dims) cout << dims << endl;
        if (print_all) cout << "Voxel Type:  ";
        if (print_all | print_type) {
            cout << imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << " " 
                 << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << endl;
        }
        if (dim3 && dims > 3) dims = 3;
        if (print_all) cout << "Size:        "; if (print_all | print_size)    { for (int i = 0; i < dims; i++) cout << imageIO->GetDimensions(i) << "\t"; cout << endl; }
        if (print_all) cout << "Spacing:     "; if (print_all | print_spacing) { for (int i = 0; i < dims; i++) cout << imageIO->GetSpacing(i) << "\t"; cout << endl; }
        if (print_all) cout << "Origin:      "; if (print_all | print_origin)  { for (int i = 0; i < dims; i++) cout << imageIO->GetOrigin(i) << "\t"; cout << endl; }
        if (print_all) cout << "Direction:   " << endl;
        if (print_all | print_direction) {
            for (int i = 0; i < dims; i++) {
                std::vector<double> dir = imageIO->GetDirection(i);
                for (int j = 0; j < dims; j++)
                    cout << dir[j] << "\t";
                cout << endl;
            }
        }
        optind++;
    }
    return EXIT_SUCCESS;
}
