/*
 *  qi_esmt.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>

#include "Util.h"
#include "ImageIO.h"
#include "IO.h"
#include "Args.h"
#include "MTFromEllipse.h"
#include "SequenceCereal.h"

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates qMT parameters from ellipse parameters.\nInputs are G, a, b.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> G_path(parser, "G_FILE", "Input G file");
    args::Positional<std::string> a_path(parser, "a_FILE", "Input a file");
    args::Positional<std::string> b_path(parser, "b_FILE", "Input b file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (in Hertz)", {'f', "f0"});
    args::ValueFlag<double> T2r_us(parser, "T2r", "T2r (in microseconds, default 12)", {"T2r"}, 12);
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag     all_residuals(parser, "RESIDUALS", "Write out all residuals", {'r',"all_resids"});
    QI::ParseArgs(parser, argc, argv, verbose);
    if (verbose) std::cout << "Opening file: " << QI::CheckPos(G_path) << std::endl;
    auto G = QI::ReadVectorImage<float>(QI::CheckPos(G_path));
    if (verbose) std::cout << "Opening file: " << QI::CheckPos(a_path) << std::endl;
    auto a = QI::ReadVectorImage<float>(QI::CheckPos(a_path));
    if (verbose) std::cout << "Opening file: " << QI::CheckPos(b_path) << std::endl;
    auto b = QI::ReadVectorImage<float>(QI::CheckPos(b_path));

    cereal::JSONInputArchive input(std::cin);
    auto seq = QI::ReadSequence<QI::SSFPMTSequence>(input, verbose);
    if (verbose) {
        std::cout << "T2r " << T2r_us.Get() << "us" << std::endl;
    }
    auto algo = std::make_shared<QI::MTFromEllipse>(seq, T2r_us.Get() * 1e-6, debug);

    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()*2);
    apply->SetOutputAllResiduals(all_residuals);
    apply->SetInput(0, G);
    apply->SetInput(1, a);
    apply->SetInput(2, b);
    if (B1) apply->SetConst(0, QI::ReadImage(B1.Get()));
    if (f0) apply->SetConst(1, QI::ReadImage(f0.Get()));
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));

    apply->SetVerbose(verbose);
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(subregion.Get()));
    }
    if (verbose) {
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;
    }
    std::string outPrefix = outarg.Get() + "EMT_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        std::string outName = outPrefix + algo->names().at(i) + QI::OutExt();
        if (verbose) std::cout << "Writing: " << outName << std::endl;
        QI::WriteImage(apply->GetOutput(i), outName);
    }
    if (verbose) std::cout << "Writing total residual." << std::endl;
    QI::WriteImage(apply->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
    if (all_residuals) {
        if (verbose) std::cout << "Writing individual residuals." << std::endl;
        QI::WriteVectorImage(apply->GetAllResidualsOutput(), outPrefix + "all_residuals" + QI::OutExt());
    }

    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
