/*
 *  SignalEquations.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  Based in part on work by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Signals/SignalEquations.h"

using namespace std;
using namespace Eigen;

namespace QI {

/******************************************************************************
 * One Component Signals
 *****************************************************************************/
VectorXcd One_MultiEcho(carrd &TE, cdbl TR, cdbl PD, cdbl T1, cdbl T2) {
	VectorXcd M = VectorXcd::Zero(TE.rows());
    M.real() = PD * (1 - exp(-TR / T1)) * exp(-TE / T2);
	return M;
}

VectorXcd One_AFI(cdbl flip, cdbl TR1, cdbl TR2, cdbl PD, cdbl T1, cdbl B1) {
	VectorXcd M = VectorXcd::Zero(2);
	const double E1 = exp(-TR1 / T1);
	const double E2 = exp(-TR2 / T1);
	const double s = sin(B1 * flip);
	const double c = cos(B1 * flip);
	M.real()[0] = PD * s * (1. - E2 + (1. - E1)*E2*c) / (1. - E1*E2*c*c);
	M.real()[1] = PD * s * (1. - E1 + (1. - E2)*E1*c) / (1. - E1*E2*c*c);
	return M;
}

} // End namespace QI
