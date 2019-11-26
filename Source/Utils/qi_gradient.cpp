/*
 *  qi_gradient.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkDerivativeImageFilter.h"

/*
 * Main
 */
int gradient_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates the gradient of an image (e.g. fieldmap).\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "FILE", "Input file");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::Log(verbose, "Reading data from: {}", QI::CheckPos(input_path));
    auto              image     = QI::ReadImage(QI::CheckPos(input_path), verbose);
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());

    QI::Log(verbose, "Calculating derivatives of image");
    auto grad = itk::DerivativeImageFilter<QI::VolumeF, QI::VolumeF>::New();
    grad->SetInput(image);
    grad->SetDirection(0);
    grad->Update();
    QI::WriteImage(grad->GetOutput(), outPrefix + "_gradx" + QI::OutExt(), verbose);
    grad->SetDirection(1);
    grad->Update();
    QI::WriteImage(grad->GetOutput(), outPrefix + "_grady" + QI::OutExt(), verbose);
    grad->SetDirection(2);
    grad->Update();
    QI::WriteImage(grad->GetOutput(), outPrefix + "_gradz" + QI::OutExt(), verbose);
    return EXIT_SUCCESS;
}
