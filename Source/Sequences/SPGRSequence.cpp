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

size_t SPGRBase::size() const {
    return FA.rows();
}

Eigen::ArrayXcd SPGRSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SPGR(p, FA, TR);
}

void SPGRSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_LOAD_DEGREES( FA );
}

void SPGRSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    QI_SEQUENCE_SAVE_DEGREES( FA );
}


/*
 * With echo-time correction
 */

Eigen::ArrayXcd SPGREchoSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SPGREcho(p, FA, TR, TE);
}

void SPGREchoSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TE", TE));
    QI_SEQUENCE_LOAD_DEGREES( FA );
}

void SPGREchoSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TE", TE));
    QI_SEQUENCE_SAVE_DEGREES( FA );
}

/*
 * With echo-time and finite-pulse corrections
 */

Eigen::ArrayXcd SPGRFiniteSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->SPGRFinite(p, FA, TR, Trf, TE);
}

void SPGRFiniteSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TE", TE));
    ar(cereal::make_nvp("Trf", Trf));
    QI_SEQUENCE_LOAD_DEGREES( FA );
}

void SPGRFiniteSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TE", TE));
    ar(cereal::make_nvp("Trf", Trf));
    QI_SEQUENCE_SAVE_DEGREES( FA );
}

} // End namespace QI