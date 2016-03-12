/*
 *  SPGR.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SIGNALS_MPRAGE_H
#define SIGNALS_MPRAGE_H

#include "QI/Signals/Common.h"

VectorXcd One_MPRAGE(cdbl flip, cdbl TR, const int Nseg, const int Nk0, carrd &TI, carrd &TD, cdbl PD, cdbl T1, cdbl B1);
Array2cd  One_MP2RAGE(const Array2d &alpha, cdbl TR, const int N, const Array3d &TD, cdbl M0, cdbl T1, cdbl B1, cdbl eta);
Array3cd  One_MP3RAGE(const Array3d &alpha, cdbl TR, const int N, const Array4d &TD, cdbl M0, cdbl T1, cdbl B1, cdbl eta);

#endif // SIGNALS_MPRAGE_H