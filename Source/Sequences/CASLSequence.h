/*
 *  CASLSequence.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_CASL_H
#define SEQUENCES_CASL_H

#include "SequenceBase.h"

namespace QI {

struct CASLSequence : QI::SequenceBase {
    double         TR, label_time;
    Eigen::ArrayXd post_label_delay;

    QI_SEQUENCE_DECLARE(CASL);
    Eigen::Index size() const override { return 1; };
};
void from_json(const json &j, CASLSequence &s);
void to_json(json &j, const CASLSequence &s);

} // End namespace QI

#endif // SEQUENCES_CASL_H