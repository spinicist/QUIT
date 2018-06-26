/*
 *  MTSat.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR sequence with saturation pulse at different offsets
 *
 */

#include "MTSatSequence.h"

namespace QI {

Eigen::Index MTSatSequence::size() const {
    return sat_f0.rows();
}

Eigen::ArrayXcd MTSatSequence::signal(const Eigen::VectorXd &p) const {
    return MT_SPGR(p, FA, TR);
}

void SPGRSequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR),
       cereal::make_nvp("FA", FA),
       cereal::make_nvp("sat_f0", sat_f0);
    QI_SEQUENCE_LOAD_DEGREES( sat_pwr );
}

void SPGRSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR),
       cereal::make_nvp("FA", FA),
       cereal::make_nvp("sat_f0", sat_f0);
    QI_SEQUENCE_SAVE_DEGREES( sat_pwr );
}


/*
 * With echo-time correction
 */

Eigen::ArrayXcd SPGREchoSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
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

Eigen::ArrayXcd SPGRFiniteSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
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