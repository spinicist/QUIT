/*
 *  qireorder.cpp
 *
 *  Created by Tobias Wood on 22/12/2015.
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>
#include <Eigen/Dense>
#include "QI/Types.h"
#include "Filters/ReorderImageFilter.h"

using namespace std;
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();

    size_t stride = 1;
    size_t blockSize = 0;
    bool verbose = false;

    const struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"stride", required_argument, 0, 's'},
        {"blocksize", required_argument, 0, 'b'},
        {"threads", required_argument, 0, 'T'},
        {0, 0, 0, 0}
    };
    const char *short_options = "hvs:b:T:h";
    const string usage {
"Usage is: reorder [options] input_file output_file \n\
\
Options:\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print more information\n\
    --stride, -s      : Set re-order stride (a good idea)\n\
    --blocksize, -b   : Only re-order within blocks\n\
    --threads, -T N   : Use N threads (default=hardware limit)\n"
};

    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
        case 'v': verbose = true; break;
        case 's': stride = atoi(optarg); break;
        case 'b': blockSize = atoi(optarg); break;
        case 'T':
            itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg));
            break;
        case 'h':
            //cout << QI::GetVersion() << endl << usage << endl;
            return EXIT_SUCCESS;
        case '?': // getopt will print an error message
            return EXIT_FAILURE;
        default:
            cout << "Unhandled option " << string(1, c) << endl;
            return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 2) {
        cout << "Incorrect number of arguments." << endl << usage << endl;
        return EXIT_FAILURE;
    }
    if (verbose) cout << "Opening input file: " << argv[optind] << endl;
    string inName(argv[optind++]);
    string outName(argv[optind++]);

    auto inFile = QI::TimeseriesReaderF::New();
    auto reorder = itk::ReorderImageFilter<QI::TimeseriesF>::New();
    auto outFile = QI::TimeseriesWriterF::New();

    inFile->SetFileName(inName);
    reorder->SetInput(inFile->GetOutput());       // Does nothing unless stride set
    if (stride > 1)
        reorder->SetStride(stride);
    if (blockSize > 0)
        reorder->SetBlockSize(blockSize);
    outFile->SetInput(reorder->GetOutput());
    outFile->SetFileName(outName);
    outFile->Update();
    if (verbose)
        cout << "Finished" << endl;
    return EXIT_SUCCESS;
}
