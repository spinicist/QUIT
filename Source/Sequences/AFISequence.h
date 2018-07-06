/*
 *  AFI.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_AFI_H
#define SEQUENCES_AFI_H

#include "SequenceBase.h"

namespace QI {

struct AFISequence : SequenceBase {
    double FA, TR1, TR2;
    Eigen::Index size() const override { return 2; }
    QI_SEQUENCE_DECLARE( AFI )
};

} // End namespace QI

#endif // SEQUENCES_AFI_H