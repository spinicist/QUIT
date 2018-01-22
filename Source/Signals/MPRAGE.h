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

#include "Common.h"

namespace QI {

Eigen::VectorXcd One_MPRAGE(cdbl flip, cdbl TR, const int Nseg, const int Nk0, cdbl TI, cdbl TD, cdbl PD, cdbl T1, cdbl B1, cdbl eta);
Eigen::Array2cd  One_MP2RAGE(const Eigen::Array2d &alpha, cdbl TR, const int N, const Eigen::Array3d &TD, cdbl M0, cdbl T1, cdbl B1, cdbl eta);
Eigen::Array3cd  One_MP3RAGE(const Eigen::Array3d &alpha, cdbl TR, const int N, const Eigen::Array4d &TD, cdbl M0, cdbl T1, cdbl B1, cdbl eta);

} // End namespace QI

#endif // SIGNALS_MPRAGE_H
