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

struct SPGRSequence : SequenceBase {
    Eigen::ArrayXd FA;
    double         TR;
    SPGRSequence(Eigen::ArrayXd const &FA, double const &TR);
    Eigen::Index size() const override;
    QI_SEQUENCE_DECLARE(SPGR);
};
void from_json(const json &j, SPGRSequence &s);
void to_json(json &j, const SPGRSequence &s);

struct SPGREchoSequence : SequenceBase {
    Eigen::ArrayXd FA;
    double         TR, TE;
    Eigen::Index   size() const override;
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
