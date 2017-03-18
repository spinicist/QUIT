/*
 *  qiesmap.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "args.hxx"

#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Args.h"
#include "QI/Algorithms/Banding.h"
#include "QI/Algorithms/Ellipse.h"
#include "QI/Algorithms/EllipseFit.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("A utility for calculating T1,T2,PD and f0 maps from SSFP data.\nInput must be a single complex image with at least 6 phase increments.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> ssfp_path(parser, "SSFP_FILE", "Input SSFP file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (h/c/f)", {'a', "algo"}, 'h');
    args::ValueFlag<int> ph_incs(parser, "INCS", "Number of phase increments (default is 6).", {"ph_incs"}, 6);
    args::Flag ph_order(parser, "PH_FIRST", "Data order is phase, then flip-angle (default opposite).", {"ph_order"});
    args::Flag alt_order(parser, "ALTERNATE", "Opposing phase-incs alternate (default is 2 blocks)", {"alternate"});
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;
    
    if (verbose) cout << "Opening file: " << QI::CheckPos(ssfp_path) << endl;
    auto data = QI::ReadVectorImage<complex<float>>(QI::CheckPos(ssfp_path));
    auto seq = make_shared<QI::SSFPEcho>(cin, prompt);
    if (verbose) cout << *seq;
    shared_ptr<QI::ESAlgo> algo;
    switch (algorithm.Get()) {
    case 'h': algo = make_shared<QI::HyperEllipse>(seq, debug, ph_order); break;
    case 'c': algo = make_shared<QI::ConstrainedEllipse>(seq, debug, ph_order, alt_order); break;
    case 'f': algo = make_shared<QI::FitEllipse>(seq, debug, ph_order, alt_order); break;
    }
    auto apply = QI::ApplyVectorXFVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    apply->SetInput(0, data);
    if (mask) {
        if (verbose) std::cout << "Reading mask: " << mask.Get() << std::endl;
        apply->SetMask(QI::ReadImage(mask.Get()));
    }
    if (B1) {
        if (verbose) std::cout << "Reading B1 map: " << B1.Get() << std::endl;
        apply->SetConst(0, QI::ReadImage(B1.Get()));
    }
    apply->SetVerbose(verbose);
    if (subregion) {
        apply->SetSubregion(QI::RegionOpt(args::get(subregion)));
    }
    if (verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    std::string outPrefix = args::get(outarg) + "ES_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        std::string outName = outPrefix + algo->names().at(i) + QI::OutExt();
        if (verbose) std::cout << "Writing: " << outName << std::endl;
        QI::WriteVectorImage(apply->GetOutput(i), outName);
    }
    
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
