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
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkTransformFileWriter.h"

#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Args.h"

using namespace std;

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::ArgParser args{argc, argv,
        "Usage is: qiaffine input [output] [transforms]\n"
        "Applies simple affine transformations to images by manipulating the header\n"
        "transforms. If an output file is not specified, the input file will be\n"
        "overwritten.",
        {{"help", 'h', "Display the help message and quit", false},
         {"verbose", 'v', "Print more information", false},
         {"center", 'c', "Set the origin to the center of the image", false},
         {"scale", 's', "Scale image by a factor of S", true},
         {"offX", '\0', "Translate origin by D in X direction", true},
         {"offY", '\0', "Translate origin by D in Y direction", true},
         {"offZ", '\0', "Translate origin by D in Z direction", true},
         {"rotX", '\0', "Rotate about X axis by N degrees", true},
         {"rotY", '\0', "Rotate about Y axis by N degrees", true},
         {"rotZ", '\0', "Rotate about Z axis by N degrees", true}}
    };

    bool verbose = args.option_present("verbose");
    bool center = args.option_present("center");
    float scale = args.option_value("scale", 1.0);
    float rotX = args.option_value("rotX", 0.0);
    float rotY = args.option_value("rotY", 0.0);
    float rotZ = args.option_value("rotZ", 0.0);
    float offX = args.option_value("offX", 0.0);
    float offY = args.option_value("offY", 0.0);
    float offZ = args.option_value("offZ", 0.0);
    std::string tfmFile = args.option_value("tfm", std::string{""});
    std::deque<const std::string> nonopts = args.nonoptions();
    if ((nonopts.size() == 0) || (nonopts.size() > 2)) {
        std::cerr << "Incorrect number of arguments, use -h to see usage." << std::endl;
        return EXIT_FAILURE;
    }

    auto image = QI::ReadImage<QI::SeriesF>(nonopts[0]);

    QI::SeriesF::DirectionType fullDir = image->GetDirection();
    QI::SeriesF::SpacingType fullSpacing = image->GetSpacing();
    QI::SeriesF::PointType fullOrigin = image->GetOrigin();
    QI::SeriesF::SizeType fullSize = image->GetLargestPossibleRegion().GetSize();
    QI::VolumeF::DirectionType direction;
    QI::VolumeF::SpacingType spacing;

    typedef itk::CenteredAffineTransform<double, 3> TAffine; 
    TAffine::OutputVectorType origin;
    QI::VolumeF::SizeType size;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            direction[i][j] = fullDir[i][j];
        }
        origin[i] = fullOrigin[i];
        spacing[i] = fullSpacing[i];
        size[i] = fullSize[i];
    }
    auto img_tfm = TAffine::New();
    itk::Versor<double> vd; vd.Set(direction);
    auto vt = itk::VersorRigid3DTransform<double>::New();
    vt->SetRotation(vd);
    img_tfm->Compose(vt);
    img_tfm->Scale(spacing);
    img_tfm->Translate(origin);

    auto tfm = TAffine::New();
    if (scale != 1.0) {
        if (verbose) cout << "Scaling by factor " << scale << endl;
        tfm->Scale(scale);
    }
    if (rotX != 0.0) {
        if (verbose) cout << "Rotating image by " << rotX << " around X axis." << endl;
        tfm->Rotate(1,2,rotX * M_PI / 180.0);
    }
    if (rotY != 0.0) {
        if (verbose) cout << "Rotating image by " << rotY << " around X axis." << endl;
        tfm->Rotate(2,0,rotY * M_PI / 180.0);
    }
    if (rotZ != 0.0) {
        if (verbose) cout << "Rotating image by " << rotZ << " around X axis." << endl;
        tfm->Rotate(0,1,rotZ * M_PI / 180.0);
    }
    itk::Versor<double>::VectorType offset;
    if (center) {
        for (int i = 0; i < 3; i++) {
            offset[i] = origin[i]-spacing[i]*size[i] / 2;
        }
    } else {
        offset[0] = offX;
        offset[1] = offY;
        offset[2] = offZ;
    }
    tfm->Translate(-offset);

    if (tfmFile != "") { // Output the transform file
        auto writer = itk::TransformFileWriterTemplate<double>::New();
        writer->SetInput(tfm);
        writer->SetFileName(tfmFile);
        writer->Update();
    }

    img_tfm->Compose(tfm);
    itk::CenteredAffineTransform<double, 3>::MatrixType fmat = img_tfm->GetMatrix();
    for (int i = 0; i < 3; i++) {
        fullOrigin[i] = img_tfm->GetOffset()[i];
    }
    for (int j = 0; j < 3; j++) {
        double scale = 0.;
        for (int i = 0; i < 3; i++) {
            scale += fmat[i][j]*fmat[i][j];
        }
        scale = sqrt(scale);
        
        fullSpacing[j] = scale;
        for (int i = 0; i < 3; i++) {
            fullDir[i][j] = fmat[i][j] / scale;
        }
    }
    image->SetDirection(fullDir);
    image->SetOrigin(fullOrigin);
    image->SetSpacing(fullSpacing);
    // Write out the edited file
    if (nonopts.size() == 2) {
        QI::WriteImage(image, nonopts[1]);
    } else {
        QI::WriteImage(image, nonopts[0]);
    }
    return EXIT_SUCCESS;
}
