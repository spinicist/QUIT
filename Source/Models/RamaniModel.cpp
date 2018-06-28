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

#include "RamaniModel.h"
#include "CerealMacro.h"

using namespace std;
using namespace Eigen;

namespace QI {
namespace Model {

Ramani::Ramani(cereal::JSONInputArchive &in) {
    QI_CLOAD(in, lineshape);
}

string Ramani::Name() const { return "qMT"; }
size_t Ramani::nParameters() const { return 9; }
const vector<string> & Ramani::ParameterNames() const {
    static vector<string> n{"PD", "T1_f", "T2_f", "T1_b", "T2_b", "k_bf", "f_b", "f0", "B1"};
    return n;
}



ArrayXXd Ramani::Bounds(const FieldStrength f) const {
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

ArrayXd Ramani::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 1.5, 0.05, 1.0, 12e-6, 1.0, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 1.5, 0.05, 1.0, 12e-6, 1.0, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool Ramani::ValidParameters(cvecd &params) const {
    // Negative values for everything except f0 make no sense
    if ((params.array().head(7) < 0.0).any()) {
        return false;
    } else {
        return true;
    }
}


Eigen::VectorXd Ramani::MTSat(cvecd &p, cdbl &FA, cdbl &TR, carrd &satf0, carrd &satflip, const RFPulse &pulse) const {
    std::cout << "START MTSat" << std::endl;
    cdbl &PD = p[0];
    cdbl &T1_f = p[1];
    cdbl &T2_f = p[2];
    cdbl &T1_b = p[3];
    cdbl &T2_b = p[4];
    cdbl &k_bf = p[5];
    cdbl &f_b  = p[6];
    cdbl &f0   = p[7];
    cdbl &B1   = p[8];

    std::cout << "Lineshape" << std::endl;
    ArrayXd W = satflip.square()*this->lineshape->value(satf0, T2_b);
    std::cout << "Calc" << std::endl;
    cdbl R1f = 1. / T1_f;
    cdbl R1r = 1. / T1_b;
    
    // F is M0r/M0b
    cdbl F = (1 - f_b) / f_b;
    cdbl kr = k_bf/F;
    carrd S = PD * F * ( R1r*kr/R1f + W + R1r + kr ) /
                    ( k_bf*(R1r + W) + (1.0 + (satflip/(2.*M_PI*satf0)).square()*(T1_f/T2_f))*(W+R1r+kr));
    std::cout << "END MTSat" << std::endl;
    return S;
}

} // End namespace Model
} // End namespace QI