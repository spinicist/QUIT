/*
 *  Direct2Algo.h
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_DIRECT2_H
#define QI_ELLIPSE_DIRECT2_H

#include <Eigen/Dense>
#include "EllipseAlgo.h"

namespace QI {

class Direct2Algo : public EllipseAlgo {
protected:
    Eigen::ArrayXd apply_internal(const Eigen::ArrayXcf &input, const double flip, const double TR, const Eigen::ArrayXd &phi, const bool debug, float &residual) const override;
public:
    Direct2Algo(std::shared_ptr<QI::SSFPEchoSequence> &seq, bool debug) : EllipseAlgo(seq, debug) {};
    size_t numOutputs() const override { return 9; }
    const std::vector<std::string> & names() const override {
        static std::vector<std::string> _names = {"G_ie", "a_ie", "b_ie", "f0_ie", 
                                                  "f_m", "a_m", "b_m", "f0_m",
                                                  "phi_rf"};
        return _names;
    }
};
    
} // End namespace QI

#endif // QI_ELLIPSE_DIRECT_H