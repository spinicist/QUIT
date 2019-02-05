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

namespace QI {

struct MPRAGESequence : SequenceBase {
    double TR, FA, eta, TI, TD;
    int    ETL, k0;
    QI_SEQUENCE_DECLARE(MPRAGE);
    Eigen::Index size() const override;
};

struct MP2RAGESequence : SequenceBase {
    double         TR;
    int            SegLength, k0;
    Eigen::Array2d FA;
    Eigen::Array3d TD;

    QI_SEQUENCE_DECLARE(MP2RAGE);
    Eigen::Index size() const override;
};

} // End namespace QI

#endif // SEQUENCES_MPRAGE_H
