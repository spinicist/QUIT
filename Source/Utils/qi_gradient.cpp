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
int gradient_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "FILE", "Input file");
    args::ValueFlag<std::string>  outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    Parse(parser);
    QI::Info("Reading data from: {}", QI::CheckPos(input_path));
    auto              image     = QI::ReadImage(QI::CheckPos(input_path));
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());

    QI::Info("Calculating derivatives of image");
    auto grad = itk::DerivativeImageFilter<QI::VolumeF, QI::VolumeF>::New();
    grad->SetInput(image);
    grad->SetDirection(0);
    grad->Update();
    QI::WriteImage(grad->GetOutput(), outPrefix + "_gradx" + QI::OutExt());
    grad->SetDirection(1);
    grad->Update();
    QI::WriteImage(grad->GetOutput(), outPrefix + "_grady" + QI::OutExt());
    grad->SetDirection(2);
    grad->Update();
    QI::WriteImage(grad->GetOutput(), outPrefix + "_gradz" + QI::OutExt());
    return EXIT_SUCCESS;
}
