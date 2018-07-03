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
#include <Eigen/Dense>
#include "JSON.h"
#include "Spline.h"

namespace QI {
namespace Lineshapes {

struct Lineshape {
    virtual std::string name() const = 0;
    virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const = 0;
    virtual rapidjson::Value jsonify(rapidjson::Document::AllocatorType &a) const = 0;
    // virtual void unserialize(rapidjson::Document &doc) const = 0;
};

struct Gaussian : Lineshape {
    virtual std::string name() const override { return "Gaussian"; };
    virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const override;
    virtual rapidjson::Value jsonify(rapidjson::Document::AllocatorType &a) const override;;
    // void save(rapidjson::Document &doc) const override;
};

struct Splineshape : Lineshape {
    double T2b_nominal = 1e-6;
    Eigen::ArrayXd frequencies, values;
    SplineInterpolator m_spline;
    Splineshape(const Eigen::ArrayXd &f0, const Eigen::ArrayXd &vals, const double T2b);
    Splineshape() = default;
    virtual std::string name() const override { return "Splineshape"; };
    virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const override;
    virtual rapidjson::Value jsonify(rapidjson::Document::AllocatorType &a) const override;
    // void save(cereal::JSONOutputArchive &ar) const override;
};

} // End namespace Lineshapes
} // End namespace QI

#endif // LINESHAPE_H
