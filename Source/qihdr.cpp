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

#include <iostream>

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Option.h"

using namespace std;

const string usage {
"Usage is: qihdr [options] input[s] \n\
\n\
By default, a summary of the header is printed. If options below are specified,\n\
only those parts of the header will be printed. Multiple files can be input,\n\
in which case the header info is written for each in order.\n\
\n"};


//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {

    QI::OptionList opts(usage);
    QI::Switch print_direction('d', "direction", "Print the image direction/orientation", opts);
    QI::Switch print_origin('o', "origin", "Print the the origin", opts);
    QI::Switch print_spacings('S',"spacings", "Print the voxel spacings/sizes", opts);
    QI::Option<int> print_spacing(0, 'P',"spacing", "Print a specific spacing/size", opts);
    QI::Switch print_sizes('s',"sizes", "Print the dimension/matrix sizes", opts);
    QI::Switch print_type('T',"dtype", "Print the data type", opts);
    QI::Switch print_dims('D',"dims","Print the number of dimensions", opts);
    QI::Switch dim3('3',"3D","Treat input as 3D (discard higher dimensions)", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() < 1) {
        std::cerr << opts << std::endl;
        std::cerr << "No input filename specified." << std::endl;
        return EXIT_FAILURE;
    }
    bool print_all = !(*print_direction || *print_origin || *print_spacings || *print_spacing ||
                       *print_sizes || *print_type || *print_dims);
    for (const string& fname : nonopts) {
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(fname.c_str(), itk::ImageIOFactory::ReadMode);
        if (!imageIO) {
            cerr << "Could not open: " << string(fname) << endl;
            break;
        }
        imageIO->SetFileName(string(fname));
        imageIO->ReadImageInformation();
        size_t dims = imageIO->GetNumberOfDimensions();
        if (print_all) cout << "File:        " << string(fname) << endl;
        if (print_all) cout << "Dimension:   "; if (print_all || *print_dims) cout << dims << endl;
        if (print_all) cout << "Voxel Type:  ";
        if (print_all | *print_type) {
            cout << imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << " " 
                 << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << endl;
        }
        if (*dim3 && dims > 3) dims = 3;
        if (print_all) cout << "Size:        "; if (print_all || *print_sizes)    { for (int i = 0; i < dims; i++) cout << imageIO->GetDimensions(i) << "\t"; cout << endl; }
        if (print_all) cout << "Spacing:     "; if (print_all || *print_spacings) { for (int i = 0; i < dims; i++) cout << imageIO->GetSpacing(i) << "\t"; cout << endl; }
        if (print_all) cout << "Origin:      "; if (print_all || *print_origin)   { for (int i = 0; i < dims; i++) cout << imageIO->GetOrigin(i) << "\t"; cout << endl; }
        if (print_all) cout << "Direction:   " << endl;
        if (print_all | *print_direction) {
            for (int i = 0; i < dims; i++) {
                std::vector<double> dir = imageIO->GetDirection(i);
                for (int j = 0; j < dims; j++)
                    cout << dir[j] << "\t";
                cout << endl;
            }
        }
        if (print_spacing.set()) {
            if (*print_spacing > -1 && *print_spacing < dims) {
                cout << imageIO->GetSpacing(*print_spacing) << endl;
            } else {
                cerr << "Invalid dimension " << *print_spacing << " for image " << fname << endl;
            }
        }
    }
    return EXIT_SUCCESS;
}
