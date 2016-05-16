/*
 *  qimask.cpp
 *
 *  Created by Tobias Wood on 2015/09/22.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <complex>
#include <getopt.h>

#include "itkComplexToModulusImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "QI/Types.h"
#include "QI/Util.h"

using namespace std;

const string usage{
"A tool to generate an image mask using an Otsu filter followed by hole-filling.\n\
This is a simple wrapper around itkOtsuThresholdImageFilter.\n\
Usage is: qimask [options] input \n\
\
For multi-volume input the first volume is used by default.\n\
Complex data can be used with -x (will take magnitude first).\n\
Output will be an image called [input]_mask.\n\
Options:\n\
    --verbose, -v     : Print more messages.\n\
    --out, -o path    : Change the output prefix.\n\
    --vol, -V N       : Use volume N. -1 will select the last volume\n\
    --complex, -x     : Input data is complex, take magnitude first\n\
    --fillh, -F N     : Fill with radius N\n\
    --threads, -T N   : Use a maximum of N threads.\n"
};

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},
    {"out",     required_argument, 0, 'o'},
    {"vol",     required_argument, 0, 'V'},
    {"complex", no_argument, 0,'x'},
    {"fillh",   required_argument, 0, 'F'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
static const char *short_options = "vo:V:xF:T:h";

int main(int argc, char **argv) {
    int indexptr = 0, c, volume = 0, fillh_radius = 0;
    string outName = "";
    bool verbose = false, is_complex = false;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
        case 'v': verbose = true; break;
        case 'o': outName = optarg; break;
        case 'V': volume = atoi(optarg); break;
        case 'x': is_complex = true; break;
        case 'F': fillh_radius = stoi(optarg); break;
        case 'T': {
            int numThreads = stoi(optarg);
            if (numThreads > 0)
                itk::MultiThreader::SetGlobalDefaultNumberOfThreads(numThreads);
        } break;
        case 'h':
            cout << QI::GetVersion() << endl << usage << endl;
            return EXIT_SUCCESS;
        case '?': // getopt will print an error message
            return EXIT_FAILURE;
        default:
            cout << "Unhandled option " << string(1, c) << endl;
            return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 1) {
        cout << QI::GetVersion() << endl << usage << endl;
        return EXIT_FAILURE;
    }
    string fname(argv[optind]);
    if (verbose) cout << "Opening input file " << fname << endl;
    if (outName == "")
        outName = fname.substr(0, fname.find(".nii")) + "_mask";
    QI::TimeseriesF::Pointer vols = ITK_NULLPTR;
    if (is_complex) {
        auto inFile = QI::TimeseriesReaderXF::New();
        inFile->SetFileName(fname);
        auto magFilter = itk::ComplexToModulusImageFilter<QI::TimeseriesXF, QI::TimeseriesF>::New();
        magFilter->SetInput(inFile->GetOutput());
        magFilter->Update();
        vols = magFilter->GetOutput();
        vols->DisconnectPipeline();
    } else {
        auto inFile = QI::TimeseriesReaderF::New();
        inFile->SetFileName(fname);
        inFile->Update();
        vols = inFile->GetOutput();
        vols->DisconnectPipeline();
    }

    // Extract one volume to process
    auto vol = itk::ExtractImageFilter<QI::TimeseriesF, QI::ImageF>::New();
    auto region = vols->GetLargestPossibleRegion();
    if (volume == -1) {
        volume = region.GetSize()[3] - 1;
    } else if (volume >= region.GetSize()[3]) {
        cerr << "Specified volume was past end of image, using first." << endl;
        volume = 0;
    }
    region.GetModifiableSize()[3] = 0;
    region.GetModifiableIndex()[3] = volume;
    vol->SetExtractionRegion(region);
    vol->SetInput(vols);
    vol->SetDirectionCollapseToSubmatrix();


    // Use an Otsu Threshold filter to generate the mask
    if (verbose) cout << "Generating Otsu mask" << endl;
    auto otsuFilter = itk::OtsuThresholdImageFilter<QI::ImageF, QI::ImageF>::New();
    otsuFilter->SetInput(vol->GetOutput());
    otsuFilter->SetOutsideValue(1);
    otsuFilter->SetInsideValue(0);
    otsuFilter->Update();
    
    QI::ImageF::Pointer finalMask = ITK_NULLPTR;
    if (fillh_radius > 0) {
        if (verbose) cout << "Filling holes" << endl;
        auto fillHoles = itk::VotingBinaryIterativeHoleFillingImageFilter<QI::ImageF>::New();
        itk::VotingBinaryIterativeHoleFillingImageFilter<QI::ImageF>::InputSizeType radius;
        radius.Fill(2);
        fillHoles->SetInput(otsuFilter->GetOutput());
        fillHoles->SetRadius(radius);
        fillHoles->SetMajorityThreshold(2); // Corresponds to (rad^3-1)/2 + 2 threshold
        fillHoles->SetBackgroundValue(0.0);
        fillHoles->SetForegroundValue(1.0);
        fillHoles->SetMaximumNumberOfIterations(3);
        fillHoles->Update();
        finalMask = fillHoles->GetOutput();
        finalMask->DisconnectPipeline();
    } else {
        finalMask = otsuFilter->GetOutput();
    }
    if (verbose) cout << "Saving mask" << endl;
    QI::WriteImage(finalMask, outName + QI::OutExt());
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

