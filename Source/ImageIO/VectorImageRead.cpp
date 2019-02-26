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

#include "ImageIO.h"
#include "ImageToVectorFilter.h"
#include "Log.h"
#include "itkImageFileReader.h"
#include <string>

namespace QI {

template <typename TVectorImg>
auto ReadImage(const std::string &path, const bool verbose) -> typename TVectorImg::Pointer {

    using TPixel    = typename TVectorImg::InternalPixelType;
    using TSeries   = itk::Image<TPixel, 4>;
    using TReader   = itk::ImageFileReader<TSeries>;
    using TToVector = itk::ImageToVectorFilter<TSeries>;

    auto file = TReader::New();
    file->SetFileName(path);
    QI::Log(verbose, "Reading image: {}", path);
    file->Update();

    auto convert = TToVector::New();
    convert->SetInput(file->GetOutput());
    QI::Log(verbose, "Converting to vector image");
    convert->Update();
    typename TVectorImg::Pointer vols = convert->GetOutput();
    if (!vols) {
        QI::Fail("Failed to read image: {}", path);
    }
    vols->DisconnectPipeline();
    return vols;
}

template auto ReadImage<QI::VectorVolumeF>(const std::string &path, const bool verbose)
    -> QI::VectorVolumeF::Pointer;
template auto ReadImage<QI::VectorVolumeXF>(const std::string &path, const bool verbose)
    -> QI::VectorVolumeXF::Pointer;

} // namespace QI

#endif // QUIT_IMAGEIO_H