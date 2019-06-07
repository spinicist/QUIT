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

namespace QI {

struct SPGRBase : SequenceBase {
    Eigen::ArrayXd FA;
    Eigen::Index   size() const override;
};

struct SPGRSequence : SPGRBase {
    double TR;
    QI_SEQUENCE_DECLARE(SPGR);
};
void from_json(const json &j, SPGRSequence &s);
void to_json(json &j, const SPGRSequence &s);

struct SPGREchoSequence : SPGRBase {
    double TR, TE;
    QI_SEQUENCE_DECLARE(SPGREcho);
};
void from_json(const json &j, SPGREchoSequence &s);
void to_json(json &j, const SPGREchoSequence &s);

// struct SPGRFiniteSequence : SPGRBase {
//     double TR, TE, Trf;
//     QI_SEQUENCE_DECLARE(SPGRFinite);
// };

} // End namespace QI

#endif // SEQUENCES_SPGR_H
