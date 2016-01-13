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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "Filters/VectorToImageFilter.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ReorderImageFilter.h"
#include "Filters/ReorderVectorFilter.h"

namespace QI {

typedef itk::Image<unsigned char, 3> ImageUC;
typedef itk::Image<int, 3> ImageI;

typedef itk::Image<float, 3> ImageF;
typedef itk::Image<float, 4> TimeseriesF;
typedef itk::VectorImage<float, 3> VectorImageF;
typedef itk::Image<std::complex<float>, 3> ImageXF;
typedef itk::Image<std::complex<float>, 4> TimeseriesXF;
typedef itk::VectorImage<std::complex<float>, 3> VectorImageXF;

typedef itk::Image<double, 3> ImageD;
typedef itk::Image<double, 4> TimeseriesD;
typedef itk::VectorImage<double, 3> VectorImageD;
typedef itk::Image<std::complex<double>, 3> ImageXD;
typedef itk::Image<std::complex<double>, 4> TimeseriesXD;
typedef itk::VectorImage<std::complex<double>, 3> VectorImageXD;

typedef itk::ImageFileReader<ImageF> ReadImageF;
typedef itk::ImageFileWriter<ImageF> WriteImageF;
typedef itk::ImageFileReader<TimeseriesF> ReadTimeseriesF;
typedef itk::ImageFileWriter<TimeseriesF> WriteTimeseriesF;
typedef itk::ReorderImageFilter<TimeseriesF> ReorderTimeseriesF;
typedef itk::ImageToVectorFilter<TimeseriesF> TimeseriesToVectorF;
typedef itk::VectorToImageFilter<VectorImageF> VectorToTimeseriesF;
typedef itk::ReorderVectorFilter<VectorImageF> ReorderVectorF;

typedef itk::ImageFileReader<ImageD> ReadImageD;
typedef itk::ImageFileWriter<ImageD> WriteImageD;

typedef itk::ImageFileReader<ImageXF> ReadImageXF;
typedef itk::ImageFileWriter<ImageXF> WriteImageXF;
typedef itk::ImageFileReader<TimeseriesXF> ReadTimeseriesXF;
typedef itk::ImageFileWriter<TimeseriesXF> WriteTimeseriesXF;
typedef itk::ReorderImageFilter<TimeseriesXF> ReorderTimeseriesXF;
typedef itk::ImageToVectorFilter<TimeseriesXF> TimeseriesToVectorXF;
typedef itk::VectorToImageFilter<VectorImageXF> VectorToTimeseriesXF;
typedef itk::ReorderVectorFilter<VectorImageXF> ReorderVectorXF;

} // End namespace QI

#endif // define QUIT_TYPES
