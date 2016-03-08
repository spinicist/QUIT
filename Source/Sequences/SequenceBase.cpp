/*
 *  SequenceBase.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Sequences/SequenceBase.h"

/******************************************************************************
 * SequenceBase
 *****************************************************************************/
ostream& operator<<(ostream& os, const SequenceBase& s) {
    s.write(os);
    return os;
}