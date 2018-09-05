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
namespace Lineshapes {

Lineshape *LineshapeFromJSON(rapidjson::Value &json) {
    rapidjson::Value::MemberIterator ls = json.MemberBegin();
    if (ls->name.GetString() == "Gaussian"s) { return new Gaussian(ls->value); }
    else if (ls->name.GetString() == "Splineshape"s) { return new Splineshape(ls->value); }
    else {
        QI_FAIL("Unknown lineshape: " << ls->name.GetString());
    }
}

// Eigen::ArrayXd Gaussian::value(const Eigen::ArrayXd &df0, const double T2b) const {
//     return sqrt(M_PI_2) * T2b * exp(-pow(2.0*M_PI*df0*T2b,2.0)/2.0);
// }

rapidjson::Value Gaussian::toJSON(rapidjson::Document::AllocatorType &/* Unused */) const {
    rapidjson::Value value(rapidjson::kObjectType);
    return value;
}

Gaussian::Gaussian(const rapidjson::Value &/* Unused */) {

}

Splineshape::Splineshape(const Eigen::ArrayXd &f0, const Eigen::ArrayXd &vals, const double T2b) {
    this->T2b_nominal = T2b;
    this->frequencies = f0;
    this->values = vals;
    this->m_spline = SplineInterpolator(f0, vals);
}

// Eigen::ArrayXd Splineshape::value(const Eigen::ArrayXd &f, const double T2b) const {
//     auto scale = T2b / T2b_nominal;
//     auto sf = f * scale;
//     return m_spline(sf) * scale;
// }

rapidjson::Value Splineshape::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value value(rapidjson::kObjectType);
    value.AddMember("T2b_nominal", T2b_nominal, a);
    value.AddMember("frequencies", ArrayToJSON(frequencies, a), a);
    value.AddMember("value", ArrayToJSON(values, a), a);
    return value;
}

Splineshape::Splineshape(const rapidjson::Value &val) {
    T2b_nominal = val["T2b_nominal"].GetDouble();
    frequencies = ArrayFromJSON(val, "frequencies");
    values = ArrayFromJSON(val, "values");
}

} // End namespace Lineshapes
} // End namespace QI
