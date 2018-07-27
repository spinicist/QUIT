/*
 *  ImageIO.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_IMAGEIO_H
#define QUIT_IMAGEIO_H

#include <string>
#include "ImageTypes.h"

namespace QI {

template<typename TImg = QI::VolumeF>
extern auto ReadImage(const std::string &path, const bool verbose = false) -> typename TImg::Pointer;

template<typename TImg = QI::VolumeF>
extern auto ReadMagnitudeImage(const std::string &path, const bool verbose = false) -> typename TImg::Pointer;

template<typename TImg>
extern void WriteImage(const TImg *ptr, const std::string &path);

template<typename TImg>
extern void WriteImage(const itk::SmartPointer<TImg> ptr, const std::string &path);

template<typename TImg>
extern void WriteMagnitudeImage(const TImg *ptr, const std::string &path);

template<typename TImg>
extern void WriteMagnitudeImage(const itk::SmartPointer<TImg> ptr, const std::string &path);

template<typename TImg>
extern void WriteScaledImage(const TImg *img, const QI::VolumeF *simg, const std::string &path);

template<typename TImg>
extern void WriteScaledImage(const itk::SmartPointer<TImg> &ptr, const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path);

template<typename TPixel = float>
extern auto ReadVectorImage(const std::string &path, const bool verbose = false) -> typename itk::VectorImage<TPixel, 3>::Pointer;

template<typename TVImg>
extern void WriteVectorImage(const TVImg *img, const std::string &path);

template<typename TVImg>
extern void WriteVectorImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path);

template<typename TVImg>
extern void WriteVectorMagnitudeImage(const TVImg *img, const std::string &path);

template<typename TVImg>
extern void WriteVectorMagnitudeImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path);

template<typename TVImg>
extern void WriteScaledVectorImage(const TVImg *img, const QI::VolumeF *simg, const std::string &path);

template<typename TVImg>
extern void WriteScaledVectorImage(const itk::SmartPointer<TVImg> &ptr, const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path);

} // End namespace QUIT

#endif // QUIT_IMAGEIO_H