/*
 *  SPGR.h
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

#ifndef SEQUENCES_SPGR_H
#define SEQUENCES_SPGR_H

#include "SequenceBase.h"
#include "cereal/cereal.hpp"
#include "EigenCereal.h"

namespace QI {

struct SPGRBase : SequenceBase {
    Eigen::ArrayXd FA;
    size_t size() const override;
};

struct SPGRSequence : SPGRBase {
    double TR;
    
    QI_SEQUENCE_DECLARE(SPGR);
};

struct SPGREchoSequence : SPGRBase {
    double TR, TE;
    Eigen::ArrayXd FA;
    QI_SEQUENCE_DECLARE(SPGREcho);
};

struct SPGRFiniteSequence : SPGRBase {
    double TR, TE, Trf;
    Eigen::ArrayXd FA;
    QI_SEQUENCE_DECLARE(SPGRFinite);
};

} // End namespace QI

#endif // SEQUENCES_SPGR_H
