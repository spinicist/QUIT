/*
 *  qi_bet.cpp
 *
 *  Created by Tobias Wood on 2015/06/03.
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include "args.hxx"
#include "ceres/ceres.h"

#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Args.h"

#include "itkImageMomentsCalculator.h"
#include "itkRegularSphereMeshSource.h"
#include "itkTriangleMeshToBinaryImageFilter.h"

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates a brain mask.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    QI::ParseArgs(parser, argc, argv);

    if (verbose) std::cout << "Reading input file: " << QI::CheckPos(input_path) << std::endl;
    auto input_image = QI::ReadImage(QI::CheckPos(input_path));

    /*
     * Calculate CoG etc.
     */
    auto moments = itk::ImageMomentsCalculator<QI::VolumeF>::New();
    moments->SetImage(input_image);
    moments->Compute();
    // ITK seems to put a negative sign on the CoG
    std::cout << moments->GetCenterOfGravity() << std::endl;
    std::cout << -moments->GetCenterOfGravity() << std::endl;
    std::cout << moments->GetCentralMoments() << std::endl;
    std::cout << moments->GetPrincipalMoments() << std::endl;
    std::cout << moments->GetPrincipalAxes() << std::endl;

    typedef itk::Mesh<double, 3> TMesh;
    typedef itk::RegularSphereMeshSource<TMesh> TEllipsoid;

    auto initial_ellipsoid = TEllipsoid::New();
    initial_ellipsoid->SetCenter(moments->GetCenterOfGravity());
    TEllipsoid::VectorType scale;
    scale.Fill(5.0);
    initial_ellipsoid->SetScale(scale);
    initial_ellipsoid->Update();
    typedef itk::TriangleMeshToBinaryImageFilter<TMesh, QI::VolumeF> TMeshToImage;
    auto mesh_to_image = TMeshToImage::New();
    mesh_to_image->SetInput(initial_ellipsoid->GetOutput());
    mesh_to_image->SetInfoImage(input_image);
    mesh_to_image->SetInsideValue(1);
    mesh_to_image->Update();
    std::string output_path = QI::StripExt(input_path.Get()) + "_mask" + QI::OutExt();
    QI::WriteImage(mesh_to_image->GetOutput(), output_path);
    return EXIT_SUCCESS;
}
