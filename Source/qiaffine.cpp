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
    QI::Switch center('c',"center","Set the origin to the center of the image", opts);
    QI::Option<std::string> tfmFile("", 't', "tfm","Save ITK transform file to specified file", opts);
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
    QI::SeriesF::SizeType fullSize = image->GetLargestPossibleRegion().GetSize();
    QI::VolumeF::DirectionType direction;
    QI::VolumeF::SpacingType spacing;
    itk::CenteredAffineTransform<double, 3>::OutputVectorType origin;
    QI::VolumeF::SizeType size;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            direction[i][j] = fullDir[i][j];
        }
        origin[i] = fullOrigin[i];
        spacing[i] = fullSpacing[i];
        size[i] = fullSize[i];
    }
    auto tfm = itk::CenteredAffineTransform<double, 3>::New();
    itk::Versor<double> vd; vd.Set(direction);
    auto vt = itk::VersorRigid3DTransform<double>::New();
    vt->SetRotation(vd);
    tfm->Compose(vt);
    tfm->Scale(spacing);
    tfm->Translate(origin);

    if (*scale != 1.0) {
        if (*verbose) cout << "Scaling by factor " << *scale << endl;
        tfm->Scale(*scale);
    }
    itk::Versor<double> rotate;
    if (*rotX != 0.0) {
        if (*verbose) cout << "Rotating image by " << *rotX << " around X axis." << endl;
        tfm->Rotate(1,2,*rotX * M_PI / 180.0);
    }
    if (*rotY != 0.0) {
        if (*verbose) cout << "Rotating image by " << *rotY << " around X axis." << endl;
        tfm->Rotate(2,0,*rotY * M_PI / 180.0);
    }
    if (*rotZ != 0.0) {
        if (*verbose) cout << "Rotating image by " << *rotZ << " around X axis." << endl;
        tfm->Rotate(0,1,*rotZ * M_PI / 180.0);
    }
    itk::Versor<double>::VectorType offset;
    if (*center) {
        for (int i = 0; i < 3; i++) {
            offset[i] = origin[i]-spacing[i]*size[i] / 2;
        }
    } else {
        offset[0] = *offX;
        offset[1] = *offY;
        offset[2] = *offZ;
    }
    tfm->Translate(offset);

    if (*tfmFile != "") { // Output the transform file
        auto writer = itk::TransformFileWriterTemplate<double>::New();
        writer->SetInput(tfm);
        writer->SetFileName(*tfmFile);
        writer->Update();
    }

    itk::CenteredAffineTransform<double, 3>::MatrixType fmat = tfm->GetMatrix();
    for (int i = 0; i < 3; i++) {
        fullOrigin[i] = tfm->GetOffset()[i];
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
