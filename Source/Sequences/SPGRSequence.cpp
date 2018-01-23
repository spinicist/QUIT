/*
 *  SPGR.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR / FLASH / FFE Sequences
 *
 */

#include "SPGRSequence.h"

namespace QI {

size_t SPGRSequence::size() const {
    return FA.rows();
}

Eigen::ArrayXcd SPGRSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SPGR(p, FA, TR);
}

/*
 * With echo-time correction
 */

size_t SPGREchoSequence::size() const {
    return FA.rows();
}

Eigen::ArrayXcd SPGREchoSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SPGREcho(p, FA, TR, TE);
}

/*
 * With echo-time and finite-pulse corrections
 */

size_t SPGRFiniteSequence::size() const {
    return FA.rows();
}

Eigen::ArrayXcd SPGRFiniteSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SPGRFinite(p, FA, TR, Trf, TE);
}

} // End namespace QI