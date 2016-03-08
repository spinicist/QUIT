/*
 *  SteadyStateSequence.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_STEADYSTATE_H
#define SEQUENCES_STEADYSTATE_H

#include "Sequences/SequenceBase.h"

class SteadyState : public SequenceBase {
    public:
        SteadyState();
        SteadyState(const ArrayXd &flip, const double TR);

        virtual size_t size() const override { return m_flip.rows(); }
};

#endif // SEQUENCES_STEADYSTATE_H