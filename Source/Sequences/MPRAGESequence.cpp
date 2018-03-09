/*
 *  MPRAGESequence.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MPRAGESequence.h"

namespace QI {

size_t MPRAGESequence::size() const { return 1; }

Eigen::ArrayXcd MPRAGESequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const {
    return m->MPRAGE(par, FA, TR, ETL, k0, eta, TI, TD);
}

void MPRAGESequence::load(cereal::JSONInputArchive &ar) {
    double FA_deg;
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TI", TI));
    ar(cereal::make_nvp("TD", TD));
    ar(cereal::make_nvp("eta", eta));
    ar(cereal::make_nvp("FA", FA_deg));
    ar(cereal::make_nvp("ETL", ETL));
    ar(cereal::make_nvp("k0", k0));
    FA = FA_deg * M_PI / 180.;
}

void MPRAGESequence::save(cereal::JSONOutputArchive &ar) const {
    double FA_deg = FA * 180. / M_PI;
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TI", TI));
    ar(cereal::make_nvp("TD", TD));
    ar(cereal::make_nvp("eta", eta));
    ar(cereal::make_nvp("FA", FA_deg));
    ar(cereal::make_nvp("ETL", ETL));
    ar(cereal::make_nvp("k0", k0));
}

/*
 * MP2RAGE
 */

size_t MP2RAGESequence::size() const { return 2; }

Eigen::ArrayXcd MP2RAGESequence::signal(const std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    QI_FAIL("Not implemented");
}

Eigen::ArrayXcd MP2RAGESequence::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP2RAGE(FA, TR, ETL, TD, M0, T1, B1, eta);
}

void MP2RAGESequence::load(cereal::JSONInputArchive &ar) {
    QI_SEQUENCE_LOAD( TR );
    QI_SEQUENCE_LOAD( TD );
    QI_SEQUENCE_LOAD( ETL );
    QI_SEQUENCE_LOAD_DEGREES( FA );
}

void MP2RAGESequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TD", TD));
    ar(cereal::make_nvp("ETL", ETL));
    QI_SEQUENCE_SAVE_DEGREES( FA );
}

/*
 * MP3RAGE
 */

size_t MP3RAGESequence::size() const { return 3; }

Eigen::ArrayXcd MP3RAGESequence::signal(const std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    QI_FAIL("Not implemented");
}

Eigen::ArrayXcd MP3RAGESequence::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP3RAGE(FA, TR, ETL, TD, M0, T1, B1, eta);
}

void MP3RAGESequence::load(cereal::JSONInputArchive &ar) {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TD", TD));
    ar(cereal::make_nvp("ETL", ETL));
    QI_SEQUENCE_LOAD_DEGREES( FA );
}

void MP3RAGESequence::save(cereal::JSONOutputArchive &ar) const {
    ar(cereal::make_nvp("TR", TR));
    ar(cereal::make_nvp("TD", TD));
    ar(cereal::make_nvp("ETL", ETL));
    QI_SEQUENCE_SAVE_DEGREES( FA );
}

} // End namespace QI
