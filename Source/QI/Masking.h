/*
 *  Masking.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_MASKING_H
#define QUIT_MASKING_H

#include <limits>

#include "QI/Types.h"

#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"

namespace QI {

VolumeI::Pointer ThresholdMask(const QI::VolumeF::Pointer &img, const float thresh);
VolumeI::Pointer OtsuMask(const QI::VolumeF::Pointer &img);

} // End namespace QI

#include "QI/Masking.hxx"

#endif