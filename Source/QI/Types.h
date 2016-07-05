/*
 * Types.h
 *
 * Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_TYPES
#define QUIT_TYPES

#include "itkImage.h"
#include "itkVectorImage.h"

#include "Filters/VectorToImageFilter.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ReorderImageFilter.h"
#include "Filters/ReorderVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"

namespace QI {

typedef itk::Image<unsigned char, 3> VolumeUC;
typedef itk::Image<int, 3> VolumeI;

typedef itk::Image<float, 3> VolumeF;
typedef itk::Image<float, 4> SeriesF;
typedef itk::VectorImage<float, 3> VectorVolumeF;
typedef itk::Image<std::complex<float>, 3> VolumeXF;
typedef itk::Image<std::complex<float>, 4> SeriesXF;
typedef itk::VectorImage<std::complex<float>, 3> VectorVolumeXF;

typedef itk::Image<double, 3> VolumeD;
typedef itk::Image<double, 4> SeriesD;
typedef itk::VectorImage<double, 3> VectorVolumeD;
typedef itk::Image<std::complex<double>, 3> VolumeXD;
typedef itk::Image<std::complex<double>, 4> SeriesXD;
typedef itk::VectorImage<std::complex<double>, 3> VectorVolumeXD;

typedef itk::ReorderImageFilter<SeriesF> ReorderSeriesF;
typedef itk::ImageToVectorFilter<SeriesF> SeriesToVectorF;
typedef itk::VectorToImageFilter<VectorVolumeF> VectorToSeriesF;
typedef itk::ReorderVectorFilter<VectorVolumeF>  ReorderVectorF;

typedef itk::ReorderImageFilter<SeriesXF> ReorderSeriesXF;
typedef itk::ImageToVectorFilter<SeriesXF> SeriesToVectorXF;
typedef itk::VectorToImageFilter<VectorVolumeXF> VectorToSeriesXF;
typedef itk::ReorderVectorFilter<VectorVolumeXF>  ReorderVectorXF;

typedef itk::ApplyAlgorithmFilter<VectorVolumeF, VolumeF, VolumeF> ApplyF;
typedef itk::ApplyAlgorithmFilter<VectorVolumeXF, VolumeF, VolumeF> ApplyXF;
typedef itk::ApplyAlgorithmFilter<VectorVolumeXF, VectorVolumeXF, VolumeF> ApplyVectorXF;
} // End namespace QI

#endif // define QUIT_TYPES
