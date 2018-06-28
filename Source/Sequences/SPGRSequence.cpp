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
#include "CerealMacro.h"
#include "CerealEigen.h"

namespace QI {

Eigen::Index SPGRBase::size() const {
    return FA.rows();
}

Eigen::ArrayXcd SPGRSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SPGR(p, FA, TR);
}

void SPGRSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD_DEGREES(ar, FA);
}

void SPGRSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE_DEGREES(ar, FA);
}


/*
 * With echo-time correction
 */

Eigen::ArrayXcd SPGREchoSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SPGREcho(p, FA, TR, TE);
}

void SPGREchoSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, TE);
    QI_CLOAD_DEGREES(ar, FA);
}

void SPGREchoSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, TE);
    QI_CSAVE_DEGREES(ar, FA);
}

/*
 * With echo-time and finite-pulse corrections
 */

Eigen::ArrayXcd SPGRFiniteSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SPGRFinite(p, FA, TR, Trf, TE);
}

void SPGRFiniteSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, TE);
    QI_CLOAD(ar, Trf);
    QI_CLOAD_DEGREES(ar, FA);
}

void SPGRFiniteSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, TE);
    QI_CSAVE(ar, Trf);
    QI_CSAVE_DEGREES(ar, FA);
}

} // End namespace QI