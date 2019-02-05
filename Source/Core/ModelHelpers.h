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

#include <cmath>

namespace QI {

/*
 *  Helper struct for converting between double/float for processing & IO
 */
template<typename DataType>
struct IOPrecision;

template<> struct IOPrecision<double> {
    using Type = float;
};

template<> struct IOPrecision<std::complex<double>> {
    using Type = std::complex<float>;
};

/*
 *  Helper struct for blocked filter output types
 */
template<bool Blocked, int ImageDimension, typename T>
struct BlockTypes {};

template<int ImageDimension, typename T>
struct BlockTypes<true, ImageDimension, T> {
    using Type = itk::VectorImage<T, ImageDimension>;
};

template<int ImageDimension, typename T>
struct BlockTypes<false, ImageDimension, T> {
    using Type = itk::Image<T, ImageDimension>;
};

template<typename T>
auto add_noise(const QI_ARRAY(T) &s, const double sigma) -> QI_ARRAY(T);

template<>
auto add_noise(const QI_ARRAY(double) &s, const double sigma) -> QI_ARRAY(double) {
    Eigen::ArrayXcd noise(s.rows());
    // Simple Box Muller transform
    Eigen::ArrayXd U = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    Eigen::ArrayXd V = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    noise.real() = (sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
    noise.imag() = (sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
    Eigen::ArrayXcd coutput(s.rows());
    coutput.real() = s + noise.real();
    coutput.imag() = noise.imag();
    return coutput.abs();
}

template<>
auto add_noise(const QI_ARRAY(std::complex<double>) &s, const double sigma) -> QI_ARRAY(std::complex<double>) {
    Eigen::ArrayXcd noise(s.rows());
    // Simple Box Muller transform
    Eigen::ArrayXd U = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    Eigen::ArrayXd V = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    noise.real() = (sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
    noise.imag() = (sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
    return s + noise;
}

} // End namespace QI

#endif // QI_MODELHELPERS_H