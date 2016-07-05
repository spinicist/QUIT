/*
 *  qMT.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Models/qMT.h"

using namespace std;
using namespace Eigen;

namespace QI {
/*
 * qMT Model (2 component)
 */

string qMT::Name() const { return "qMT"; }
size_t qMT::nParameters() const { return 9; }
const vector<string> & qMT::ParameterNames() const {
    static vector<string> n{"PD", "T1f", "T2f", "T1r", "T2r", "kf", "F", "f0", "B1"};
    return n;
}

ArrayXXd qMT::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.700, 0.030, 0.010, 0.000001, 0.5, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 2.000, 0.200, 5.000, 0.001000, 4.5, 0.999, 0.0, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.700, 0.030, 0.010, 0.000001, 0.5, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 2.000, 0.200, 5.000, 0.001000, 4.5, 0.999, 0.0, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd qMT::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 1.5, 0.05, 1.0, 12e-6, 1.0, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 1.5, 0.05, 1.0, 12e-6, 1.0, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool qMT::ValidParameters(cvecd &params) const {
    // Negative values for everything except f0 make no sense
    if ((params.array().head(7) < 0.0).any()) {
        return false;
    } else {
        return true;
    }
}

void qMT::setLineshape(TLineshape &l) { m_lineshape = l; }

VectorXcd qMT::SPGR_MT(cvecd &p, carrd &satflip, carrd &satf0, cdbl flip, cdbl TR, cdbl Trf) const {
    return scale(MT_SPGR(satflip, satf0, TR, Trf, m_lineshape, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]));
}

} // End namespace QI