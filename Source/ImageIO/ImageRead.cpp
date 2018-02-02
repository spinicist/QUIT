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
auto ReadImage(const std::string &path) -> typename TImg::Pointer {
    typedef itk::ImageFileReader<TImg> TReader;
    typename TReader::Pointer file = TReader::New();
    file->SetFileName(path);
    file->Update();
    typename TImg::Pointer img = file->GetOutput();
    if (!img) {
        QI_EXCEPTION("Failed to read file: " << path);
    }
    img->DisconnectPipeline();
    return img;
}

template<typename TImg>
auto ReadMagnitudeImage(const std::string &path) -> typename TImg::Pointer {
    typedef itk::Image<std::complex<typename TImg::PixelType>, TImg::ImageDimension> TComplex;
    auto x_img = QI::ReadImage<TComplex>(path);
    auto magFilter = itk::ComplexToModulusImageFilter<TComplex, TImg>::New();
    magFilter->SetInput(x_img);
    magFilter->Update();
    auto img = magFilter->GetOutput();
    img->DisconnectPipeline();
    return img;
}

template auto ReadImage<VolumeF>(const std::string &path) -> typename VolumeF::Pointer;
template auto ReadImage<VolumeD>(const std::string &path) -> typename VolumeD::Pointer;
template auto ReadImage<VolumeXF>(const std::string &path) -> typename VolumeXF::Pointer;
template auto ReadImage<VolumeXD>(const std::string &path) -> typename VolumeXD::Pointer;
template auto ReadImage<VolumeI>(const std::string &path) -> typename VolumeI::Pointer;
template auto ReadImage<VolumeUC>(const std::string &path) -> typename VolumeUC::Pointer;
template auto ReadImage<SeriesF>(const std::string &path) -> typename SeriesF::Pointer;
template auto ReadImage<SeriesD>(const std::string &path) -> typename SeriesD::Pointer;
template auto ReadImage<SeriesXF>(const std::string &path) -> typename SeriesXF::Pointer;
template auto ReadImage<SeriesXD>(const std::string &path) -> typename SeriesXD::Pointer;
template auto ReadMagnitudeImage<VolumeF>(const std::string &path) -> typename VolumeF::Pointer;
template auto ReadMagnitudeImage<SeriesF>(const std::string &path) -> typename SeriesF::Pointer;

} // End namespace QUIT

#endif // QUIT_IMAGEIO_H