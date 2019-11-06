/*
 *  ModelHelpers.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_MODELHELPERS_H
#define QI_MODELHELPERS_H

#include "Macro.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include <cmath>
#include <complex>

namespace QI {

/*
 *  Helper struct for converting between double/float for processing & IO
 */
template <typename DataType> struct IOPrecision;

template <> struct IOPrecision<double> { using Type = float; };

template <> struct IOPrecision<std::complex<double>> { using Type = std::complex<float>; };

/*
 *  Helper struct for blocked filter output types
 */
template <bool Blocked, int ImageDimension, typename T> struct BlockTypes {};

template <int ImageDimension, typename T> struct BlockTypes<true, ImageDimension, T> {
    using Type = itk::VectorImage<T, ImageDimension>;
};

template <int ImageDimension, typename T> struct BlockTypes<false, ImageDimension, T> {
    using Type = itk::Image<T, ImageDimension>;
};

Eigen::ArrayXd  add_noise(Eigen::ArrayXd const &s, double const sigma);
Eigen::ArrayXcd add_noise(Eigen::ArrayXcd const &s, double const sigma);

} // End namespace QI

#endif // QI_MODELHELPERS_H