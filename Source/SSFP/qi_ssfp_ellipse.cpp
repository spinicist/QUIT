/*
 *  qi_ellipse.cpp
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
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

#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "DirectAlgo.h"
#include "HyperAlgo.h"

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates the ellipse parameters G,a,b,f0 & psi0 from SSFP data.\nInput must be a single complex image with at least 6 phase increments.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> ssfp_path(parser, "SSFP_FILE", "Input SSFP file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (h)yper/(d)irect, default d", {'a', "algo"}, 'd');
    QI::ParseArgs(parser, argc, argv, verbose);
    QI_LOG(verbose, "Opening file: " << QI::CheckPos(ssfp_path));
    auto data = QI::ReadVectorImage<std::complex<float>>(QI::CheckPos(ssfp_path));
    rapidjson::Document input = QI::ReadJSON(std::cin);
    QI::SSFPEllipseSequence seq(input["SSFPEllipse"]);
    std::shared_ptr<QI::EllipseAlgo> algo;
    switch (algorithm.Get()) {
    case 'h': algo = std::make_shared<QI::HyperAlgo>(seq, debug); break;
    case 'd': algo = std::make_shared<QI::DirectAlgo>(seq, debug); break;
    }
    QI::ApplyVectorXFVectorF::Pointer apply = QI::ApplyVectorXFVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()*2);
    apply->SetInput(0, data);
    if (mask) {
        QI_LOG(verbose, "Reading mask: " << mask.Get());
        apply->SetMask(QI::ReadImage(mask.Get()));
    }
    if (B1) apply->SetConst(0, QI::ReadImage(B1.Get()));
    apply->SetVerbose(verbose);
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    }
    QI_LOG(verbose, "Processing");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    QI_LOG(verbose, "Elapsed time was " << apply->GetTotalTime() << "s" <<
                    "Writing results files.");
    std::string outPrefix;
    outPrefix = outarg.Get() + "ES_";
    for (size_t i = 0; i < algo->numOutputs(); i++) {
        std::string outName = outPrefix + algo->names().at(i) + QI::OutExt();
        QI_LOG(verbose, "Writing: " << outName );
        QI::WriteVectorImage(apply->GetOutput(i), outName);
    }
    QI_LOG(verbose, "Writing total residuals." );
    QI::WriteVectorImage(apply->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
    QI_LOG(verbose, "Finished." );
    return EXIT_SUCCESS;
}
