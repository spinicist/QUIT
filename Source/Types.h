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

namespace QI {

typedef itk::Image<float, 3> FloatImage;
typedef itk::VectorImage<float, 3> FloatVectorImage;
typedef itk::Image<float, 4> FloatTimeseries;

typedef itk::Image<std::complex<float>, 3> XFloatImage;
typedef itk::VectorImage<std::complex<float>, 3> XFloatVectorImage;
typedef itk::Image<std::complex<float>, 4> XFloatTimeseries;

typedef itk::ImageFileReader<FloatImage> ReadFloatImage;
typedef itk::ImageFileReader<FloatTimeseries> ReadFloatTimeseries;
typedef itk::ImageToVectorFilter<FloatTimeseries> FloatTimeseriesToVector;

typedef itk::ImageFileReader<XFloatImage> ReadXFloatImage;
typedef itk::ImageFileReader<XFloatTimeseries> ReadXFloatTimeseries;
typedef itk::ImageToVectorFilter<XFloatTimeseries> XFloatTimeseriesToVector;

} // End namespace QI
