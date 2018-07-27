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
#include "ImageToVectorFilter.h"
#include "ImageIO.h"
#include "Macro.h"

namespace QI {

template<typename TPixel>
auto ReadVectorImage(const std::string &path, const bool verbose) -> typename itk::VectorImage<TPixel, 3>::Pointer {
    typedef itk::Image<TPixel, 4> TSeries;
    typedef itk::VectorImage<TPixel, 3> TVector;
    typedef itk::ImageToVectorFilter<TSeries> TToVector;
    
    auto img = ReadImage<TSeries>(path);
    auto convert = TToVector::New();
    convert->SetInput(img);
    QI_LOG( verbose, "Reading image: " << path );
    convert->Update();
    typename TVector::Pointer vols = convert->GetOutput();
    vols->DisconnectPipeline();
    return vols;
}

template auto ReadVectorImage<float>(const std::string &path, const bool verbose) -> typename itk::VectorImage<float, 3>::Pointer;
template auto ReadVectorImage<std::complex<float>>(const std::string &path, const bool verbose) -> typename itk::VectorImage<std::complex<float>, 3>::Pointer;

} // End namespace QUIT

#endif // QUIT_IMAGEIO_H