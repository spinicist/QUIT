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
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageMomentsCalculator.h"

#include "Util.h"
#include "IO.h"
#include "Args.h"

using namespace std;

/*
 * Declare arguments here so they are available to pipeline
 */
args::ArgumentParser parser("Applies simple affine transformations to images by manipulating the header\n"
    "transforms. If an output file is not specified, the input file will be\n"
    "overwritten.\n"
    "http://github.com/spinicist/QUIT");

args::Positional<std::string> source_path(parser, "SOURCE", "Source file");
args::Positional<std::string> dest_path(parser, "DEST", "Destination file");

args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
args::ValueFlag<std::string> center(parser, "CENTER", "Set the origin to geometric center (geo) or (cog)", {'c', "center"});
args::ValueFlag<std::string> tfm_path(parser, "TFM", "Write out the transformation to a file", {'t', "tfm"});
args::ValueFlag<double> scale(parser, "SCALE", "Scale by a constant", {'s', "scale"});
args::ValueFlag<double> offX(parser, "OFF_X", "Translate origin in X direction", {"offX"});
args::ValueFlag<double> offY(parser, "OFF_Y", "Translate origin in Y direction", {"offY"});
args::ValueFlag<double> offZ(parser, "OFF_Z", "Translate origin in Z direction", {"offZ"});
args::ValueFlag<double> rotX(parser, "ROT_X", "Rotate about X-axis by angle (degrees)", {"rotX"});
args::ValueFlag<double> rotY(parser, "ROT_Y", "Rotate about Y-axis by angle (degrees)", {"rotY"});
args::ValueFlag<double> rotZ(parser, "ROT_Z", "Rotate about Z-axis by angle (degrees)", {"rotZ"});
args::ValueFlag<std::string> permute(parser, "PERMUTE", "Permute axes, e.g. 2,0,1. Negative values mean flip as well", {"permute"});
args::ValueFlag<std::string> flip(parser, "FLIP", "Flip an axis, e.g. 0,1,0. Occurs AFTER any permutation.", {"flip"});

template<typename TImage>
int Pipeline() {
    auto image = QI::ReadImage<TImage>(QI::CheckPos(source_path));

    // Permute if required
    if (permute) {
        auto permute_filter = itk::PermuteAxesImageFilter<TImage>::New();
        itk::FixedArray<unsigned int, TImage::ImageDimension> permute_order;
        std::istringstream iss(permute.Get());
        std::string el;
        for (int i = 0; i < 3; i++) {
            std::getline(iss, el, ',');
            permute_order[i] = std::stoi(el);
        }
        for (int i = 3; i < TImage::ImageDimension; i++) {
            permute_order[i] = i;
        }
        if (!iss)
            QI_FAIL("Failed to read permutation order: " << permute.Get());

        permute_filter->SetInput(image);
        permute_filter->SetOrder(permute_order);
        permute_filter->Update();
        image = permute_filter->GetOutput();
        image->DisconnectPipeline();
    }

    // Flip if required
    if (flip) {
        auto flip_filter = itk::FlipImageFilter<TImage>::New();
        itk::FixedArray<bool, TImage::ImageDimension> flip_axes; // Save this
        std::istringstream iss(flip.Get());
        std::string el;
        for (int i = 0; i < 3; i++) {
            std::getline(iss, el, ',');
            flip_axes[i] = (std::stoi(el) > 0);
        }
        for (int i = 3; i < TImage::ImageDimension; i++) {
            flip_axes[i] = false;
        }
        if (!iss)
            QI_FAIL("Failed to read flip: " << flip.Get());
        if (verbose) std::cout << "Flipping: " << flip_axes << std::endl;
        flip_filter->SetInput(image);
        flip_filter->SetFlipAxes(flip_axes);
        flip_filter->SetFlipAboutOrigin(false);
        flip_filter->Update();
        image = flip_filter->GetOutput();
        image->DisconnectPipeline();
    }


    typename TImage::DirectionType fullDir = image->GetDirection();
    typename TImage::SpacingType fullSpacing = image->GetSpacing();
    typename TImage::PointType fullOrigin = image->GetOrigin();
    typename TImage::SizeType fullSize = image->GetLargestPossibleRegion().GetSize();
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
    img_tfm->SetMatrix(direction);
    img_tfm->Scale(spacing);
    img_tfm->Translate(origin);
    if (verbose) std::cout << "Start transform:\n" << img_tfm->GetMatrix() << std::endl;
    auto tfm = TAffine::New();
    if (scale) {
        if (verbose) cout << "Scaling by factor " << scale.Get() << endl;
        tfm->Scale(scale.Get());
    }
    if (rotX != 0.0) {
        if (verbose) cout << "Rotating image by " << rotX.Get() << " around X axis." << endl;
        tfm->Rotate(1,2,rotX.Get() * M_PI / 180.0);
    }
    if (rotY != 0.0) {
        if (verbose) cout << "Rotating image by " << rotY.Get() << " around X axis." << endl;
        tfm->Rotate(2,0,rotY.Get() * M_PI / 180.0);
    }
    if (rotZ != 0.0) {
        if (verbose) cout << "Rotating image by " << rotZ.Get() << " around X axis." << endl;
        tfm->Rotate(0,1,rotZ.Get() * M_PI / 180.0);
    }
    itk::Versor<double>::VectorType offset; offset.Fill(0);
    if (center) {
        if (center.Get() == "geo") {
            if (verbose) std::cout << "Setting geometric center" << std::endl;
            for (int i = 0; i < 3; i++) {
                offset[i] = origin[i]-spacing[i]*size[i] / 2;
            }
        } else if (center.Get() == "cog") {
            if (verbose) std::cout << "Setting center to center of gravity" << std::endl;
            auto moments = itk::ImageMomentsCalculator<TImage>::New();
            moments->SetImage(image);
            moments->Compute();
            // ITK seems to put a negative sign on the CoG
            for (int i = 0; i < 3; i++) {
                offset[i] = moments->GetCenterOfGravity()[i];
            }
        }
        std::cout << "Translation will be: " << offset << std::endl;
        tfm->Translate(-offset);
    } else if (offX || offY || offZ) {
        offset[0] = offX.Get();
        offset[1] = offY.Get();
        offset[2] = offZ.Get();
        if (verbose) std::cout << "Translating by: " << offset << std::endl;
        tfm->Translate(-offset);
    }

    if (tfm_path) { // Output the transform file
        auto writer = itk::TransformFileWriterTemplate<double>::New();
        writer->SetInput(tfm);
        writer->SetFileName(tfm_path.Get());
        writer->Update();
    }

    img_tfm->Compose(tfm);
    itk::CenteredAffineTransform<double, 3>::MatrixType fmat = img_tfm->GetMatrix();
    if (verbose) std::cout << "Final transform:\n" << fmat << std::endl;
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
        if (verbose) std::cout << "Writing: " << dest_path.Get() << std::endl;
        QI::WriteImage(image, dest_path.Get());
    } else {
        if (verbose) std::cout << "Writing: " << source_path.Get() << std::endl;
        QI::WriteImage(image, source_path.Get());
    }
    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    QI::ParseArgs(parser, argc, argv);
    if (verbose) std::cout << "Reading header for: " << QI::CheckPos(source_path) << std::endl;
    auto header = itk::ImageIOFactory::CreateImageIO(QI::CheckPos(source_path).c_str(), itk::ImageIOFactory::ReadMode);
    if (!header)  {
        QI_EXCEPTION("Failed to read header from: " << source_path.Get());
    }
    header->SetFileName(source_path.Get());
    header->ReadImageInformation();
    auto dims  = header->GetNumberOfDimensions();
    auto dtype = header->GetComponentType();
    if (verbose) std::cout << "Datatype is " << header->GetComponentTypeAsString( dtype ) << std::endl;

    #define DIM_SWITCH( TYPE ) \
        switch (dims) { \
            case 3: return Pipeline<itk::Image< TYPE , 3 >>(); \
            case 4: return Pipeline<itk::Image< TYPE , 4 >>(); \
            default: QI_FAIL("Unsupported dimension: " << dims); return EXIT_FAILURE; \
        }

    switch (dtype) {
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE: QI_FAIL("Unknown component type in image " << source_path);
        case itk::ImageIOBase::UCHAR:  DIM_SWITCH( unsigned char ); break;
        case itk::ImageIOBase::CHAR:   DIM_SWITCH( char ); break;
        case itk::ImageIOBase::USHORT: DIM_SWITCH( unsigned short ); break;
        case itk::ImageIOBase::SHORT:  DIM_SWITCH( short ); break;
        case itk::ImageIOBase::UINT:   DIM_SWITCH( unsigned int ); break;
        case itk::ImageIOBase::INT:    DIM_SWITCH( int ); break;
        case itk::ImageIOBase::ULONG:  DIM_SWITCH( unsigned long ); break;
        case itk::ImageIOBase::LONG:   DIM_SWITCH( long ); break;
        case itk::ImageIOBase::LONGLONG: DIM_SWITCH( long long ); break;
        case itk::ImageIOBase::ULONGLONG: DIM_SWITCH( unsigned long long ); break;
        case itk::ImageIOBase::FLOAT:  DIM_SWITCH( float ); break;
        case itk::ImageIOBase::DOUBLE: DIM_SWITCH( double ); break;
    }
}
