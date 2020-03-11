/*
 *  qiaffine.cpp
 *
 *  Copyright (c) 2020 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Args.h"
#include "Util.h"
#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkTransformFileReader.h"

int affine_angle_main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Calculates the composite rotation angle of the Z-axis by a set of transforms\n"
        "http://github.com/spinicist/QUIT");

    args::PositionalList<std::string> tfm_paths(parser, "TRANSFORM", "List of transform files");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     inverse(parser, "INVERSE", "Apply inverse transforms", {'i', "inverse"});
    QI::ParseArgs(parser, argc, argv, verbose);

    using Tfm      = itk::AffineTransform<double, 3>;
    auto composite = Tfm::New();
    composite->SetIdentity();
    // using Tfm = itk::CompositeTransform<double, 3>;

    for (auto const &tfm_path : tfm_paths.Get()) {

        QI::Info(verbose, "Transform file: {}", tfm_path);

        auto reader = itk::TransformFileReader::New();
        reader->SetFileName(tfm_path);
        reader->Update();

        auto const &tfm  = *(reader->GetTransformList()->begin());
        auto const &atfm = static_cast<Tfm *>(tfm.GetPointer());
        QI::Log(verbose, "{}", *atfm);
        if (inverse) {
            auto itfm = Tfm::New();
            atfm->GetInverse(itfm);
            composite->Compose(itfm);
        } else {
            composite->Compose(atfm);
        }
    }

    QI::Log(verbose, "Composite:\n{}", composite);
    itk::CovariantVector<double> Z0;
    Z0[0] = 0.;
    Z0[1] = 0.;
    Z0[2] = 1.;

    itk::CovariantVector Z1 = composite->TransformCovariantVector(Z0);
    QI::Log(verbose, "Z0: {}\nZ1: {}", Z0, Z1);
    auto const norm = Z1.Normalize();
    QI::Log(verbose, "Z1 norm: {}", norm);
    auto const angle = acos(Z0 * Z1);
    fmt::print("{}\n", angle * 180. / M_PI);
    return EXIT_SUCCESS;
}
