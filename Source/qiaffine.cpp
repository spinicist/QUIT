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

    args::ArgumentParser parser("Applies simple affine transformations to images by manipulating the header\n"
        "transforms. If an output file is not specified, the input file will be\n"
        "overwritten.\n"
        "http://github.com/spinicist/QUIT");

    args::Positional<std::string> source_path(parser, "SOURCE", "Source file");
    args::Positional<std::string> dest_path(parser, "DEST", "Destination file");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     center(parser, "CENTER", "Set the origin to the center of the image", {'c', "center"});
    args::ValueFlag<std::string> tfm_path(parser, "TFM", "Write out the transformation to a file", {'t', "tfm"});
    args::ValueFlag<double> scale(parser, "SCALE", "Scale by a constant", {'s', "scale"});
    args::ValueFlag<double> offX(parser, "OFF_X", "Translate origin in X direction", {"offX"});
    args::ValueFlag<double> offY(parser, "OFF_Y", "Translate origin in Y direction", {"offY"});
    args::ValueFlag<double> offZ(parser, "OFF_Z", "Translate origin in Z direction", {"offZ"});
    args::ValueFlag<double> rotX(parser, "ROT_X", "Rotate about X-axis by angle (degrees)", {"rotX"});
    args::ValueFlag<double> rotY(parser, "ROT_Y", "Rotate about Y-axis by angle (degrees)", {"rotY"});
    args::ValueFlag<double> rotZ(parser, "ROT_Z", "Rotate about Z-axis by angle (degrees)", {"rotZ"});

    QI::ParseArgs(parser, argc, argv);

    auto image = QI::ReadImage<QI::SeriesF>(QI::CheckPos(source_path));

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

    if (tfm_path) { // Output the transform file
        auto writer = itk::TransformFileWriterTemplate<double>::New();
        writer->SetInput(tfm);
        writer->SetFileName(tfm_path.Get());
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
    if (dest_path) {
        QI::WriteImage(image, dest_path.Get());
    } else {
        QI::WriteImage(image, source_path.Get());
    }
    return EXIT_SUCCESS;
}
