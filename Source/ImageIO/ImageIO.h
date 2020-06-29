#pragma once
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

#include "ImageTypes.h"
#include <string>

namespace QI {

template <typename TImg = QI::VolumeF>
extern auto ReadImage(const std::string &path, const bool verbose) -> typename TImg::Pointer;

template <typename TImg = QI::VolumeF>
extern auto ReadMagnitudeImage(const std::string &path, const bool verbose) ->
    typename TImg::Pointer;

template <typename TImg>
extern void WriteImage(const TImg *ptr, const std::string &path, const bool verbose);

template <typename TImg>
extern void
WriteImage(const itk::SmartPointer<TImg> &ptr, const std::string &path, const bool verbose);

template <typename TImg>
extern void WriteMagnitudeImage(const TImg *ptr, const std::string &path, const bool verbose);

template <typename TImg>
extern void WriteMagnitudeImage(const itk::SmartPointer<TImg> &ptr,
                                const std::string &            path,
                                const bool                     verbose);

template <typename TImg>
extern void WriteScaledImage(const TImg *       img,
                             const QI::VolumeF *simg,
                             const std::string &path,
                             const bool         verbose);

template <typename TImg>
extern void WriteScaledImage(const itk::SmartPointer<TImg> &       ptr,
                             const itk::SmartPointer<QI::VolumeF> &sptr,
                             const std::string &                   path,
                             const bool                            verbose);

} // namespace QI
