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
#include <cereal/archives/json.hpp>

namespace QI {
namespace Lineshapes {

struct Lineshape {
    virtual std::string name() const = 0;
    virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const = 0;
    virtual void load(cereal::JSONInputArchive &ar) = 0;
    virtual void save(cereal::JSONOutputArchive &ar) const = 0;
};

struct Gaussian : Lineshape {
    virtual std::string name() const override { return "Gaussian"; };
    virtual Eigen::ArrayXd value(const Eigen::ArrayXd &f, const double T2b) const override;
    void load(cereal::JSONInputArchive &ar) override;
    void save(cereal::JSONOutputArchive &ar) const override;
};

} // End namespace Lineshapes
} // End namespace QI

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::Lineshapes::Lineshape> const &ls);
void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::Lineshapes::Lineshape> &ls);

} // End namespace cereal

#endif // LINESHAPE_H
