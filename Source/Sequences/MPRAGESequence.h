/*
 *  MPRAGESequence.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_MPRAGE_H
#define SEQUENCES_MPRAGE_H

#include "SequenceBase.h"
#include "Macro.h"

namespace QI {

struct MPRAGESequence : SequenceBase {
    double TR = 0.0, FA = 0.0, eta = 0.0, TI = 0.0, TD = 0.0;
    int ETL = 0, k0 = 0;
    QI_SEQUENCE_DECLARE(MPRAGE);
    Eigen::Index size() const override;
};

struct MP2RAGESequence : SequenceBase {
    double TR = 0.0;
    int ETL = 0;
    Eigen::ArrayXd FA;
    Eigen::Array3d TD = Eigen::Array3d::Zero();

    QI_SEQUENCE_DECLARE(MP2RAGE);
    Eigen::Index size() const override;
    Eigen::ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
};

} // End namespace QI

#endif // SEQUENCES_MPRAGE_H
