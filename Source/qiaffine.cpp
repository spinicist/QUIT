/*
 *  qiaffine.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>

#include "itkImage.h"
#include "itkVersor.h"
#include "itkChangeInformationImageFilter.h"

#include "QI/Util.h"
#include "QI/Types.h"

using namespace std;

int main(int argc, char **argv) {
    const string usage {
    "Usage is: qiaffine input [output] [transforms] \n\
    \n\
    Applies simple affine transformations to images by manipulating the header\n\
    transforms. If an output file is not specified, the input file will be\n\
    overwritten.\n\
    \n\
    Transform Options:\n\
        --scale S : Scale all axes by a factor of S\n\
        --rotX N  : Rotate about the X axis by N degrees\n\
        --rotY N  : Rotate about the Y axis by N degrees\n\
        --rotZ N  : Rotate about the Z axis by N degrees\n\
        --origin 'X Y Z' : Set the origin to X Y Z\n\
    \n\
    Other Options:\n\
        --help, -h    : Print this message\n\
        --verbose, -v : Print more messages\n\
    \n"};

    int verbose = false;
    const struct option long_options[] =
    {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"scale", required_argument, 0, 'S'},
        {"rotX", required_argument, 0, 'X'},
        {"rotY", required_argument, 0, 'Y'},
        {"rotZ", required_argument, 0, 'Z'},
        {"origin", required_argument, 0, 'O'},
        {0, 0, 0, 0}
    };
    const char* short_options = "hv";
    int c, index_ptr = 0;

    // Do one pass to get options in the right order
    while ((c = getopt_long(argc, argv, short_options, long_options, &index_ptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'h':
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case 'S': case 'X': case 'Y': case 'Z': case 'O':
            case 0: break; // A flag
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
            cout << "Unhandled option " << string(1, c) << endl;
            return EXIT_FAILURE;
        }
    }

    if (((argc - optind) == 0) || ((argc - optind) > 2)) {
        cout << QI::GetVersion() << endl << usage << endl;
        cout << "Incorrect number of arguments" << endl;
        return EXIT_FAILURE;
    }

    // Now read in the input image
    auto reader = itk::ImageFileReader<QI::SeriesF>::New();
    reader->SetFileName(argv[optind]);
    reader->Update();
    auto image = reader->GetOutput();
    auto writer = itk::ImageFileWriter<QI::SeriesF>::New();
    if ((argc - optind) == 2) {
        optind++;
    }
    writer->SetFileName(argv[optind]);

    // Now reset optind and actually process
    QI::SeriesF::DirectionType fullDir = image->GetDirection();
    QI::SeriesF::SpacingType fullSpacing = image->GetSpacing();
    QI::SeriesF::PointType fullOrigin = image->GetOrigin();
    QI::VolumeF::DirectionType direction;
    QI::VolumeF::SpacingType spacing;
    QI::VolumeF::PointType origin;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            direction[i][j] = fullDir[i][j];
        }
        origin[i] = fullOrigin[i];
        spacing[i] = fullSpacing[i];
    }
    optind = 1;
    while ((c = getopt_long(argc, argv, short_options, long_options, &index_ptr)) != -1) {
        switch (c) {
            // Already handled these
            case 'v': case 'h': case 0: break;
            case 'S': {
                const double fac = atof(optarg);
                if (verbose) cout << "Scaling by factor " << fac << endl;
                spacing = spacing * fac;
            } break;
            case 'X': {
                if (verbose) cout << "Rotating image by " << string(optarg) << " around X axis." << endl;
                const double radians = atof(optarg) * vnl_math::pi / 180.0;
                itk::Versor<double> rotate; rotate.SetRotationAroundX(radians);
                direction = rotate.GetMatrix() * direction;
                origin = rotate.GetMatrix() * origin;
            } break;
            case 'Y': {
                if (verbose) cout << "Rotating image by " << string(optarg) << " around Y axis." << endl;
                const double radians = atof(optarg) * vnl_math::pi / 180.0;
                itk::Versor<double> rotate; rotate.SetRotationAroundY(radians);
                direction = rotate.GetMatrix() * direction;
                origin = rotate.GetMatrix() * origin;
            } break;
            case 'Z': {
                if (verbose) cout << "Rotating image by " << string(optarg) << " around Z axis." << endl;
                const double radians = atof(optarg) * vnl_math::pi / 180.0;
                itk::Versor<double> rotate; rotate.SetRotationAroundZ(radians);
                direction = rotate.GetMatrix() * direction;
                origin = rotate.GetMatrix() * origin;
            } break;
            case 'O': {
                string option(optarg);
                if (verbose) cout << "Setting origin to " << option << endl;
                stringstream optstream(option);
                optstream >> origin;
            } break;
        }
    }
    auto changeInfo = itk::ChangeInformationImageFilter<QI::SeriesF>::New();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            fullDir[i][j] = direction[i][j];
        }
        fullOrigin[i] = origin[i];
        fullSpacing[i] = spacing[i];
    }
    changeInfo->SetOutputDirection(fullDir);
    changeInfo->SetOutputOrigin(fullOrigin);
    changeInfo->SetOutputSpacing(fullSpacing);
    changeInfo->ChangeDirectionOn();
    changeInfo->ChangeOriginOn();
    changeInfo->ChangeSpacingOn();
    changeInfo->SetInput(image);
    writer->SetInput(changeInfo->GetOutput());
    writer->Update();
    return EXIT_SUCCESS;
}
