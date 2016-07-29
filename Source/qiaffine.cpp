/*
 *  qiaffine.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "itkImage.h"
#include "itkVersor.h"
#include "itkVersorTransform.h"
#include "itkTransformFileWriter.h"
#include "itkChangeInformationImageFilter.h"

#include "QI/Util.h"
#include "QI/Types.h"
#include "QI/Option.h"

using namespace std;

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts("Usage is: qiaffine input [output] [transforms]\n\nApplies simple affine transformations to images by manipulating the header\ntransforms. If an output file is not specified, the input file will be\noverwritten.");
    QI::Option<float> scale(1,'\0',"scale","Scale axes by a factor of S", opts);
    QI::Option<float> rotX(0,'\0',"rotX","Rotate about X axis by N degrees", opts);
    QI::Option<float> rotY(0,'\0',"rotY","Rotate about Y axis by N degrees", opts);
    QI::Option<float> rotZ(0,'\0',"rotZ","Rotate about Z axis by N degrees", opts);
    QI::Option<float> offX(0,'\0',"offX","Offset X", opts);
    QI::Option<float> offY(0,'\0',"offY","Offset Y", opts);
    QI::Option<float> offZ(0,'\0',"offZ","Offset Z", opts);
    QI::Option<std::string> tfmFile("", 't', "tfm","Save ITK transform file", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if ((nonopts.size() == 0) || (nonopts.size() > 2)) {
        std::cerr << opts << std::endl;
        std::cerr << "Incorrect number of arguments" << std::endl;
        return EXIT_FAILURE;
    }

    // Now read in the input image
    auto reader = itk::ImageFileReader<QI::SeriesF>::New();
    reader->SetFileName(nonopts[0]);
    reader->Update();
    auto image = reader->GetOutput();

    QI::SeriesF::DirectionType fullDir = image->GetDirection();
    QI::SeriesF::SpacingType fullSpacing = image->GetSpacing();
    QI::SeriesF::PointType fullOrigin = image->GetOrigin();
    QI::VolumeF::DirectionType direction;
    QI::VolumeF::SpacingType spacing;
    QI::VolumeF::PointType origin;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            direction[i][j] = fullDir[i][j];
        }
        origin[i] = fullOrigin[i];
        spacing[i] = fullSpacing[i];
    }
    
    if (*scale != 1.0) {
        if (*verbose) cout << "Scaling by factor " << *scale << endl;
        spacing = spacing * (*scale);
    }
    itk::Versor<double> rotate;
    if (*rotX != 0.0) {
        if (*verbose) cout << "Rotating image by " << string(optarg) << " around X axis." << endl;
        itk::Versor<double> temp; temp.SetRotationAroundX(*rotX * M_PI / 180.0);
        rotate *= temp;
    }
    if (*rotY != 0.0) {
        if (*verbose) cout << "Rotating image by " << string(optarg) << " around X axis." << endl;
        itk::Versor<double> temp; temp.SetRotationAroundY(*rotY * M_PI / 180.0);
        rotate *= temp;
    }
    if (*rotZ != 0.0) {
        if (*verbose) cout << "Rotating image by " << string(optarg) << " around X axis." << endl;
        itk::Versor<double> temp; temp.SetRotationAroundZ(*rotZ * M_PI / 180.0);
        rotate *= temp;
    }
    itk::Versor<double>::VectorType offset;
    offset[0] = *offX;
    offset[1] = *offY;
    offset[2] = *offZ;
    direction = rotate.GetMatrix() * direction;
    origin = rotate.GetMatrix() * origin;
    
    /*
        string option(optarg);
        if (verbose) cout << "Setting origin to " << option << endl;
        stringstream optstream(option);
        optstream >> origin;
    */

    if (*tfmFile != "") { // Only output the transform file
        auto tfm = itk::VersorTransform<double>::New();
        tfm->SetCenter(origin);
        tfm->SetRotation(rotate);
        tfm->SetOffset(offset);
        auto writer = itk::TransformFileWriterTemplate<double>::New();
        writer->SetInput(tfm);
        writer->SetFileName(*tfmFile);
        writer->Update();
    } else { // Write out the edited file
        auto changeInfo = itk::ChangeInformationImageFilter<QI::SeriesF>::New();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                fullDir[i][j] = direction[i][j];
            }
            fullOrigin[i] = origin[i];
            fullSpacing[i] = spacing[i];
        }
        changeInfo->SetOutputDirection(fullDir);
        changeInfo->SetOutputOrigin(fullOrigin);
        changeInfo->SetOutputSpacing(fullSpacing);
        changeInfo->ChangeDirectionOn();
        changeInfo->ChangeOriginOn();
        changeInfo->ChangeSpacingOn();
        changeInfo->SetInput(image);
        if (nonopts.size() == 2) {
            QI::WriteImage(changeInfo->GetOutput(), nonopts[1]);
        } else {
            QI::WriteImage(changeInfo->GetOutput(), nonopts[0]);
        }
    }
    return EXIT_SUCCESS;
}
