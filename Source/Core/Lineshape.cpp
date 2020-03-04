/*
 *  Lineshape.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Lineshape.h"
#include "Macro.h"

using namespace std::string_literals;

namespace QI {

InterpLineshape::InterpLineshape(const double          fmin,
                                 const double          fstep,
                                 const int             fcount,
                                 const Eigen::ArrayXd &vals,
                                 const double          T2) :
    T2_nominal{T2},
    freq_min{fmin}, freq_step{fstep}, freq_count{fcount}, values{vals} {

    grid         = std::make_shared<ceres::Grid1D<double>>(&values[0], 0, values.rows());
    interpolator = std::make_shared<ceres::CubicInterpolator<ceres::Grid1D<double>>>(*grid);
}

} // End namespace QI

namespace nlohmann {

QI::InterpLineshape adl_serializer<QI::InterpLineshape>::from_json(const json &j) {
    auto T2nom = j.at("T2_nominal").get<double>();
    auto fmin  = j.at("freq_min").get<double>();

    auto fstep  = j.at("freq_step").get<double>();
    auto fcount = j.at("freq_count").get<int>();
    auto vals   = QI::ArrayFromJSON<double>(j, "values");

    return QI::InterpLineshape(fmin, fstep, fcount, vals, T2nom);
}

void adl_serializer<QI::InterpLineshape>::to_json(json &j, QI::InterpLineshape l) {
    j = json{{"T2_nominal", l.T2_nominal},
             {"freq_min", l.freq_min},
             {"freq_step", l.freq_step},
             {"freq_count", l.freq_count},
             {"values", l.values}};
}

} // namespace nlohmann
