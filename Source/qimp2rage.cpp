/*
 *  qimp2rage.cpp
 *
 *  Created by Tobias Wood on 2015/08/24.
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

#include "Types.h"
#include "Util.h"

#include "itkBinaryFunctorImageFilter.h"
#include "itkMaskImageFilter.h"

using namespace std;

const string usage{
"Usage is: qimp2rage [options] input \n\
\
Input must contain two volumes and be complex-valued\n\
Options:\n\
    --verbose, -v   : Print more messages\n\
    --mask, -m file : Mask input with specified file.\n\
    --out, -o path  : Add a prefix to the output filenames.\n"
};

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},
    {"mask",    required_argument, 0, 'm'},
    {"out",     required_argument, 0, 'o'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
static const char *short_options = "vm:o:T:h";

template<class T> class MP2RAGE {
public:
    MP2RAGE() {}
    ~MP2RAGE() {}
    bool operator!=(const MP2RAGE &) const { return false; }
    bool operator==(const MP2RAGE &other) const { return !(*this != other); }

    inline T operator()(const complex<T> &ti1,
                              const complex<T> &ti2) const
    {
        const T a1 = abs(ti1);
        const T a2 = abs(ti2);
        return real((conj(ti1)*ti2)/(a1*a1 + a2*a2));
    }
};


int main(int argc, char **argv) {
    int indexptr = 0, c;
    string outPrefix = "";
    bool verbose = false;
    QI::ReadImageF::Pointer maskFile = ITK_NULLPTR;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'm':
                cout << "Reading mask." << endl;
                maskFile = QI::ReadImageF::New();
                maskFile->SetFileName(optarg);
                break;
            case 'o':
                outPrefix = optarg;
                cout << "Output prefix will be: " << outPrefix << endl;
                break;
            case 'T': itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg)); break;
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
    if ((argc - optind) != 1) {
        cout << usage << endl;
        return EXIT_FAILURE;
    }
    string fname(argv[optind]);
    if (verbose) cout << "Opening input file " << fname << endl;
    if (outPrefix == "")
        outPrefix = fname.substr(0, fname.find(".nii"));
    auto inFile = QI::ReadTimeseriesXF::New();
    inFile->SetFileName(fname);
    inFile->Update();
    vector<itk::ExtractImageFilter<QI::TimeseriesXF, QI::ImageXF>::Pointer> vols(3);

    auto region = inFile->GetOutput()->GetLargestPossibleRegion();
    region.GetModifiableSize()[3] = 0;

    for (int i = 0; i < 3; i++) {
        region.GetModifiableIndex()[3] = i;
        vols[i] = itk::ExtractImageFilter<QI::TimeseriesXF, QI::ImageXF>::New();
        vols[i]->SetExtractionRegion(region);
        vols[i]->SetInput(inFile->GetOutput());
        vols[i]->SetDirectionCollapseToSubmatrix();
    }

    for (int i1 = 0; i1 < 3; i1++) {
        for (int i2 = (i1 + 1); i2 < 3; i2++) {
            auto mp2rage_filter = itk::BinaryFunctorImageFilter<QI::ImageXF, QI::ImageXF, QI::ImageF, MP2RAGE<float>>::New();
            mp2rage_filter->SetInput1(vols[i1]->GetOutput());
            mp2rage_filter->SetInput2(vols[i2]->GetOutput());

            auto mask_filter = itk::MaskImageFilter<QI::ImageF, QI::ImageF, QI::ImageF>::New();
            if (maskFile) {
                mask_filter->SetInput1(mp2rage_filter->GetOutput());
                mask_filter->SetMaskImage(maskFile->GetOutput());
                mask_filter->Update();
                QI::writeResult<QI::ImageF>(mask_filter->GetOutput(), outPrefix + "_C" + to_string(i1) + to_string(i2) + QI::OutExt());
            } else {
                mp2rage_filter->Update();
                QI::writeResult<QI::ImageF>(mp2rage_filter->GetOutput(), outPrefix + "_C" + to_string(i1) + to_string(i2) + QI::OutExt());
            }
        }
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

