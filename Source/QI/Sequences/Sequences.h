/*
 *  Sequences.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SEQUENCES_H
#define SEQUENCES_SEQUENCES_H

#include "QI/Sequences/SteadyStateSequence.h"
#include "QI/Sequences/MPRAGESequence.h"
#include "QI/Sequences/SpinEcho.h"

namespace QI {
    std::shared_ptr<SequenceBase> ReadSequence(std::istream& istr, const bool prompt);
} // End namespace QI

#endif // SEQUENCES_SEQUENCES_H