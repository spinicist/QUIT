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

#include "Types.h"

namespace QI {

VolumeI::Pointer ThresholdMask(const QI::VolumeF::Pointer &img,
                               const float lower,
                               const float upper = std::numeric_limits<float>::infinity());
VolumeI::Pointer OtsuMask(const QI::VolumeF::Pointer &img);
std::vector<float> FindLabels(const QI::VolumeI::Pointer &mask, const int size_threshold, const int to_keep, QI::VolumeI::Pointer &labels);

} // End namespace QI

#endif