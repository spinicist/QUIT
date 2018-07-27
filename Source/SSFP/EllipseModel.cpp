/*
 *  EllipseModel.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "EllipseModel.h"

using namespace std::literals;

namespace QI {

std::array<const std::string, EllipseModel::NV> EllipseModel::varying_names{{"G"s, "a"s, "b"s, "theta_0"s, "phi_rf"s}};
std::array<const std::string, EllipseModel::NF> EllipseModel::fixed_names{};
const QI_ARRAYN(double, EllipseModel::NF) EllipseModel::fixed_defaults{};

auto EllipseModel::signal(const QI_ARRAYN(double, NV) &v,
                          const QI_ARRAYN(double, NF) &/*Unused*/,
                          const QI::SSFPSequence *s) const -> QI_ARRAY(std::complex<double>)
{
    const double &G = v[0];
    const double &a = v[1];
    const double &b = v[2];

    if (b > 2.*a/(1. + a*a)) {
        return QI_ARRAY(std::complex<double>)::Zero(s->PhaseInc.rows());
    }

    const double &theta0 = v[3];
    const double &psi0 = v[4];
    const auto theta = theta0 - s->PhaseInc;
    const double psi = theta0/2.0 + psi0;
    const auto cos_th = cos(theta);
    const auto sin_th = sin(theta);
    const auto cos_psi = cos(psi);
    const auto sin_psi = sin(psi);
    QI_ARRAY(std::complex<double>) result(s->PhaseInc.rows());
    result.real() = G * (cos_psi - a*(cos_th*cos_psi - sin_th*sin_psi)) / (1.0 - b*cos_th);
    result.imag() = G * (sin_psi - a*(cos_th*sin_psi + sin_th*cos_psi)) / (1.0 - b*cos_th);
    // QI_DB( G )
    // QI_DB( a )
    // QI_DB( b )
    // QI_DB( theta0 )
    // QI_DB( psi0 )
    // QI_DBVEC( s->PhaseInc )
    // QI_DBVEC( theta )
    // QI_DBVEC( cos_th )
    // QI_DBVEC( sin_th )
    // QI_DB( cos_psi )
    // QI_DB( sin_psi )
    // QI_DBVEC( result )
    return result;
}
} // End namespace QI