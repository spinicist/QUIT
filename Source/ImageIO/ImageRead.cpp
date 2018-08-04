/*
 *  ImageRead.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_IMAGEIO_H

#include <string>

#include "itkImageFileReader.h"
#include "itkComplexToModulusImageFilter.h"
#include "ImageIO.h"
#include "Macro.h"

namespace QI {

template<typename TImg>
auto ReadImage(const std::string &path, const bool verbose) -> typename TImg::Pointer {
    typedef itk::ImageFileReader<TImg> TReader;
    typename TReader::Pointer file = TReader::New();
    file->SetFileName(path);
    QI_LOG( verbose, "Reading image: " << path );
    file->Update();
    typename TImg::Pointer img = file->GetOutput();
    if (!img) {
        QI_EXCEPTION("Failed to read file: " << path);
    }
    img->DisconnectPipeline();
    return img;
}

template<typename TImg>
auto ReadMagnitudeImage(const std::string &path, const bool verbose) -> typename TImg::Pointer {
    typedef itk::Image<std::complex<typename TImg::PixelType>, TImg::ImageDimension> TComplex;
    QI_LOG( verbose, "Reading image: " << path );
    auto x_img = QI::ReadImage<TComplex>(path, verbose);
    auto magFilter = itk::ComplexToModulusImageFilter<TComplex, TImg>::New();
    magFilter->SetInput(x_img);
    QI_LOG( verbose, "Converting to magnitude" );
    magFilter->Update();
    auto img = magFilter->GetOutput();
    img->DisconnectPipeline();
    return img;
}

template auto ReadImage<VolumeF>(const std::string &path, const bool verbose) -> typename VolumeF::Pointer;
template auto ReadImage<VolumeD>(const std::string &path, const bool verbose) -> typename VolumeD::Pointer;
template auto ReadImage<VolumeXF>(const std::string &path, const bool verbose) -> typename VolumeXF::Pointer;
template auto ReadImage<VolumeXD>(const std::string &path, const bool verbose) -> typename VolumeXD::Pointer;
template auto ReadImage<VolumeI>(const std::string &path, const bool verbose) -> typename VolumeI::Pointer;
template auto ReadImage<VolumeUC>(const std::string &path, const bool verbose) -> typename VolumeUC::Pointer;
template auto ReadImage<SeriesF>(const std::string &path, const bool verbose) -> typename SeriesF::Pointer;
template auto ReadImage<SeriesD>(const std::string &path, const bool verbose) -> typename SeriesD::Pointer;
template auto ReadImage<SeriesXF>(const std::string &path, const bool verbose) -> typename SeriesXF::Pointer;
template auto ReadImage<SeriesXD>(const std::string &path, const bool verbose) -> typename SeriesXD::Pointer;
template auto ReadMagnitudeImage<VolumeF>(const std::string &path, const bool verbose) -> typename VolumeF::Pointer;
template auto ReadMagnitudeImage<SeriesF>(const std::string &path, const bool verbose) -> typename SeriesF::Pointer;

} // End namespace QUIT

#endif // QUIT_IMAGEIO_H