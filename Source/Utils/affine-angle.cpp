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

int affine_angle_main(args::Subparser &parser) {
    args::PositionalList<std::string> tfm_paths(parser, "TRANSFORM", "List of transform files");
    parser.Parse();

    using Tfm      = itk::AffineTransform<double, 3>;
    auto composite = Tfm::New();
    composite->SetIdentity();

    for (auto const &tfm_path : tfm_paths.Get()) {
        bool const        inverse = (tfm_path[0] == '^');
        std::string const path    = inverse ? tfm_path.substr(1) : tfm_path;
        QI::Info(verbose, "{}Transform file: {}", inverse ? "Inverse " : "", path);

        auto reader = itk::TransformFileReader::New();
        reader->SetFileName(path);
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
