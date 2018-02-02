/*
 *  ApplyAlgorithmFilter.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ApplyAlgorithmFilter.hxx"
#include "ImageTypes.h"

namespace itk {

template class ApplyAlgorithmFilter<QI::VectorVolumeF, QI::VectorVolumeF, QI::VolumeF, QI::VolumeF>;
template class ApplyAlgorithmFilter<QI::VectorVolumeXF, QI::VectorVolumeXF, QI::VolumeF, QI::VolumeF>;
template class ApplyAlgorithmFilter<QI::VectorVolumeXF, QI::VectorVolumeF, QI::VolumeF, QI::VolumeF>;

}