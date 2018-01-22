/*
 *  SSFP.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SSFP.h"
#include "IO.h"

namespace QI {

Eigen::ArrayXd SSFP::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFP::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFP(p, FA, TR, PhaseInc);
}

Eigen::ArrayXcd SSFPEcho::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFPEcho(p, FA, TR, PhaseInc);
}

Eigen::ArrayXd SSFPEcho::signal_magnitude(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFPEchoMagnitude(p, FA, TR, PhaseInc);
}

Eigen::ArrayXd SSFPFinite::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFPFinite::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFPFinite(p, FA, TR, Trf, PhaseInc);
}

Eigen::ArrayXd SSFP_GS::weights(const double f0) const {
    return Eigen::ArrayXd::Ones(size()); // Weight SPGR images higher than SSFP
}

Eigen::ArrayXcd SSFP_GS::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SSFP_GS(p, FA, TR);
}

} // End namespace QI
