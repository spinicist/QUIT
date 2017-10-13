/*
 *  DirectAlgo.h
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_DIRECT_H
#define QI_ELLIPSE_DIRECT_H

#include <Eigen/Dense>
#include "QI/Ellipse/EllipseAlgo.h"

namespace QI {

class DirectAlgo : public EllipseAlgo {
protected:
    Eigen::ArrayXd apply_internal(const Eigen::ArrayXcf &input, const double flip, const double TR, const Eigen::ArrayXd &phi, const bool debug, float &residual) const override;
public:
    DirectAlgo(std::shared_ptr<QI::SSFPEcho> &seq, bool debug) : EllipseAlgo(seq, debug) {};
    size_t numOutputs() const override { return 5; }
    const std::vector<std::string> & names() const override {
        static std::vector<std::string> _names = {"G", "a", "b", "f0", "phi_rf"};
        return _names;
    }
};

} // End namespace QI

#endif // QI_ELLIPSE_DIRECT_H