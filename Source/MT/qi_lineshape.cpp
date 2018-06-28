/*
 *  qi_lineshape.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood, Samuel Hurley, Erika Raven
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "Util.h"
#include "Args.h"
#include "Lineshape.h"
#include "CerealEigen.h"

int main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser("Saves lineshapes as JSON objects.\n"
                                "http://github.com/spinicist/QUIT");

    args::Positional<std::string> output_path(parser, "OUTPUT", "Output JSON file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag                    verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string>  lineshape(parser, "LINESHAPE", "Choose lineshape (Gauss/Lorentzian/Super-lorentzian)", {'l', "lineshape"});

    QI::ParseArgs(parser, argc, argv, verbose);
    if (verbose) std::cout << "Reading input from: " << QI::CheckPos(input_path) << std::endl;
    auto input = QI::ReadImage(QI::CheckPos(input_path));
    auto fit = itk::PolynomialFitImageFilter::New();
    QI::Polynomial<3> poly(order.Get());
    fit->SetInput(input);
    fit->SetPolynomial(poly);
    fit->SetRobust(robust);
    Eigen::Array3d center = Eigen::Array3d::Zero();
    if (mask_path) {
        if (verbose) std::cout << "Reading mask from: " << mask_path.Get() << std::endl;
        auto mask_image = QI::ReadImage(mask_path.Get());
        fit->SetMask(mask_image);
        auto moments = itk::ImageMomentsCalculator<QI::VolumeF>::New();
        moments->SetImage(mask_image);
        moments->Compute();
        // ITK seems to put a minus sign on CoG
        if (verbose) std::cout << "Mask CoG is: " << -moments->GetCenterOfGravity() << std::endl;
        center = QI::Eigenify(-moments->GetCenterOfGravity());
    }
    fit->SetCenter(center);
    QI::VolumeF::SizeType size = input->GetLargestPossibleRegion().GetSize();
    QI::VolumeF::SpacingType spacing = input->GetSpacing();
    for (int i = 0; i < 3; i++) spacing[i] *= size[i];
    double scale = (spacing/2).GetNorm();
    fit->SetScale(scale);
    fit->Update();
    {
        cereal::JSONOutputArchive output(std::cout);
        QI::WriteCereal(output, "center", center);
        QI::WriteCereal(output, "scale", scale);
        QI::WriteCereal(output, "coeffs", fit->GetPolynomial().coeffs());
        if (print_terms)
            QI::WriteCereal(output, "terms", fit->GetPolynomial().get_terms());
        QI::WriteCereal(output, "residual",  fit->GetResidual());
    }
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
