/*
 *  SpinEcho.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SPINECHO_H
#define SEQUENCES_SPINECHO_H

#include "SequenceBase.h"

namespace QI {

struct MultiEchoSequence : SequenceBase {
    double TR, TE1, ESP;
    int ETL;
    Eigen::ArrayXd TE;
    QI_SEQUENCE_DECLARE(MultiEcho);
    size_t size() const override;
};

} // End namespace QI

#endif // SEQUENCES_SPINECHO_H
