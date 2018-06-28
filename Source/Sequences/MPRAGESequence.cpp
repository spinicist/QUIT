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
#include "CerealMacro.h"

namespace QI {

Eigen::Index MPRAGESequence::size() const { return 1; }

Eigen::ArrayXcd MPRAGESequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &par) const {
    return m->MPRAGE(par, FA, TR, ETL, k0, eta, TI, TD);
}

void MPRAGESequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, TI);
    QI_CLOAD(ar, TD);
    QI_CLOAD(ar, eta);
    QI_CLOAD_DEGREES(ar, FA);
    QI_CLOAD(ar, ETL);
    QI_CLOAD(ar, k0);
}

void MPRAGESequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, TI);
    QI_CSAVE(ar, TD);
    QI_CSAVE(ar, eta);
    QI_CSAVE_DEGREES(ar, FA);
    QI_CSAVE(ar, ETL);
    QI_CSAVE(ar, k0);
}

/*
 * MP2RAGE
 */

Eigen::Index MP2RAGESequence::size() const { return 2; }

Eigen::ArrayXcd MP2RAGESequence::signal(const std::shared_ptr<Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not implemented");
}

Eigen::ArrayXcd MP2RAGESequence::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP2RAGE(FA, TR, ETL, TD, M0, T1, B1, eta);
}

void MP2RAGESequence::load(cereal::JSONInputArchive &ar) {
    double SegTR;
    Eigen::Array2d TI;
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, SegTR);
    QI_CLOAD(ar, TI);
    QI_CLOAD(ar, ETL);
    QI_CLOAD_DEGREES(ar, FA);
    TD[0] = TI[0];
    TD[1] = TI[1] - (ETL * TR) - TI[0];
    TD[2] = SegTR - (ETL * TR) - TI[1];
}

void MP2RAGESequence::save(cereal::JSONOutputArchive &ar) const {
    Eigen::Array2d TI{TD[0], TD[1] + (ETL * TR) + TD[0]};
    double SegTR = TI[1] + (ETL * TR) + TD[2];
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, SegTR);
    QI_CSAVE(ar, TI);
    QI_CSAVE(ar, ETL);
    QI_CSAVE_DEGREES(ar, FA);
}

/*
 * MP3RAGE
 */

Eigen::Index MP3RAGESequence::size() const { return 3; }

Eigen::ArrayXcd MP3RAGESequence::signal(const std::shared_ptr<Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not implemented");
}

Eigen::ArrayXcd MP3RAGESequence::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP3RAGE(FA, TR, ETL, TD, M0, T1, B1, eta);
}

void MP3RAGESequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar,  TD);
    QI_CLOAD(ar, ETL);
    QI_CLOAD_DEGREES(ar, FA);
}

void MP3RAGESequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar,  TD);
    QI_CSAVE(ar, ETL);
    QI_CSAVE_DEGREES(ar, FA);
}

} // End namespace QI
