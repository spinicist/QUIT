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

InterpLineshape::InterpLineshape(const double fmin, const double fstep, const int fcount,
                                 const Eigen::ArrayXd &vals, const double T2) :
    T2_nominal{T2},
    freq_min{fmin}, freq_step{fstep}, freq_count{fcount},
    values{vals},
    grid(&values[0], 0, values.rows()),
    interpolator{grid}
{}

InterpLineshape::InterpLineshape(const rapidjson::Value &val) :
    T2_nominal{QI::GetMember(val, "T2_nominal").GetDouble()},
    freq_min{QI::GetMember(val, "freq_min").GetDouble()},
    freq_step{QI::GetMember(val, "freq_step").GetDouble()},
    freq_count{QI::GetMember(val, "freq_count").GetInt()},
    values{ArrayFromJSON(val, "values")},
    grid(&values[0], 0, values.rows()),
    interpolator{grid}
{}

rapidjson::Value InterpLineshape::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value value(rapidjson::kObjectType);
    value.AddMember("T2_nominal", T2_nominal, a);
    value.AddMember("freq_min", freq_min, a);
    value.AddMember("freq_step", freq_step, a);
    value.AddMember("freq_count", freq_count, a);
    value.AddMember("values", ArrayToJSON(values, a), a);
    return value;
}

} // End namespace QI
