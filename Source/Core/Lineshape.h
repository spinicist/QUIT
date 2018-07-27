/*
 *  Lineshape.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef LINESHAPE_H
#define LINESHAPE_H

#include <functional>
#include <string>
#include <Eigen/Core>
#include "JSON.h"
#include "Spline.h"
#include "Macro.h"

namespace QI {
namespace Lineshapes {

struct Lineshape {
    virtual std::string name() const = 0;
    // virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const = 0;
    virtual rapidjson::Value toJSON(rapidjson::Document::AllocatorType &a) const = 0;
};
Lineshape *LineshapeFromJSON(rapidjson::Value &);

struct Gaussian : Lineshape {
    virtual std::string name() const override { return "Gaussian"; };
    template<typename T>
    QI_ARRAY(T) value(const Eigen::ArrayXd &f0, const T T2b) const {
        return sqrt(M_PI_2) * T2b * exp(-pow(2.0*M_PI*f0*T2b,2.0)/2.0);
    }
    virtual rapidjson::Value toJSON(rapidjson::Document::AllocatorType &a) const override;
    Gaussian() = default;
    Gaussian(const rapidjson::Value &);
};

struct Splineshape : Lineshape {
    double T2b_nominal = 1e-6;
    Eigen::ArrayXd frequencies, values;
    SplineInterpolator m_spline;
    virtual std::string name() const override { return "Splineshape"; };
    // virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const override;
    virtual rapidjson::Value toJSON(rapidjson::Document::AllocatorType &a) const override;
    Splineshape() = default;
    Splineshape(const rapidjson::Value &);
    Splineshape(const Eigen::ArrayXd &f0, const Eigen::ArrayXd &vals, const double T2b);
};

} // End namespace Lineshapes
} // End namespace QI

#endif // LINESHAPE_H
