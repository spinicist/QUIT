/*
 * ApplyTypes.h
 *
 * Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_APPLY_TYPES
#define QUIT_APPLY_TYPES

#include "ImageTypes.h"
#include "ApplyAlgorithmFilter.h"

namespace QI {

typedef itk::ApplyAlgorithmFilter<VectorVolumeF, VolumeF, VolumeF, VolumeF> ApplyF;
typedef itk::ApplyAlgorithmFilter<VectorVolumeXF, VolumeF, VolumeF, VolumeF> ApplyXF;
typedef itk::ApplyAlgorithmFilter<VectorVolumeF, VectorVolumeF, VolumeF, VolumeF> ApplyVectorF;
typedef itk::ApplyAlgorithmFilter<VectorVolumeXF, VectorVolumeXF, VolumeF, VolumeF> ApplyVectorXF;
typedef itk::ApplyAlgorithmFilter<VectorVolumeXF, VectorVolumeF, VolumeF, VolumeF> ApplyVectorXFVectorF;

} // End namespace QI

#endif // define QUIT_APPLY_TYPES
