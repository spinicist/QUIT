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

#include <iostream>
#include <Eigen/Dense>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "ImageTypes.h"
#include "Args.h"
#include "ReorderImageFilter.h"

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Reorders image series\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT", "Path to input series file");
    args::Positional<std::string> output_path(parser, "OUTPUT", "Path to output file");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<int> stride(parser, "STRIDE", "Re-order stride (will do nothing if not set)", {'s', "stride"}, 0);
    args::ValueFlag<int> blocksize(parser, "BLOCKSIZE", "Only re-order within blocks, size N (default one block)", {'b', "blocksize"}, 0);
    QI::ParseArgs(parser, argc, argv);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());
    if (verbose) std::cout << "Opening input file: " << QI::CheckPos(input_path) << std::endl;
    std::string inName(QI::CheckPos(input_path));
    std::string outName(QI::CheckPos(output_path));

    auto inFile = itk::ImageFileReader<QI::SeriesF>::New();
    auto reorder = itk::ReorderImageFilter<QI::SeriesF>::New();
    auto outFile = itk::ImageFileWriter<QI::SeriesF>::New();

    inFile->SetFileName(inName);
    reorder->SetInput(inFile->GetOutput());       // Does nothing unless stride set
    if (stride) reorder->SetStride(stride.Get());
    if (blocksize) reorder->SetBlockSize(blocksize.Get());
    outFile->SetInput(reorder->GetOutput());
    outFile->SetFileName(outName);
    outFile->Update();
    if (verbose) std::cout << "Finished" << std::endl;
    return EXIT_SUCCESS;
}
