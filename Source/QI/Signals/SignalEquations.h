/*
 *  SignalEquations.h
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SIGNAL_EQUATIONS_H
#define SIGNAL_EQUATIONS_H

#include "QI/Signals/SPGR.h"
#include "QI/Signals/SSFP.h"
#include "QI/Signals/MPRAGE.h"

VectorXcd One_MultiEcho(carrd &TE, cdbl TR, cdbl PD, cdbl T1, cdbl T2);
VectorXcd One_AFI(cdbl flip, cdbl TR1, cdbl TR2, cdbl PD, cdbl T1, cdbl B1);

#endif
