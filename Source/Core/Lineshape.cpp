/*
 *  Lineshape.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Lineshape.h"
#include "CerealMacro.h"
#include "CerealEigen.h"

namespace QI {
namespace Lineshapes {

Eigen::ArrayXd Gaussian::value(const Eigen::ArrayXd &df0, const double T2b) const {
    return sqrt(M_PI_2) * T2b * exp(-pow(2.0*M_PI*df0*T2b,2.0)/2.0);
}

void Gaussian::load(cereal::JSONInputArchive &/* Unused */) {
}

void Gaussian::save(cereal::JSONOutputArchive &/* Unused */) const {
}

Splineshape::Splineshape(const Eigen::ArrayXd &f0, const Eigen::ArrayXd &vals, const double T2b) {
    this->T2b_nominal = T2b;
    this->frequencies = f0;
    this->values = vals;
    this->m_spline = SplineInterpolator(f0, vals);
}

Eigen::ArrayXd Splineshape::value(const Eigen::ArrayXd &f, const double T2b) const {
    auto scale = T2b / T2b_nominal;
    auto sf = f * scale;
    return m_spline(sf) * scale;
}

void Splineshape::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, T2b_nominal);
    QI_CLOAD(ar, frequencies);
    QI_CLOAD(ar, values);
}

void Splineshape::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, T2b_nominal);
    QI_CSAVE(ar, frequencies);
    QI_CSAVE(ar, values);
}

} // End namespace Lineshapes
} // End namespace QI

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::Lineshapes::Lineshape> const &ls) {
    #define QI_SAVE( NAME ) \
        (ls->name() == #NAME ) { ar(cereal::make_nvp(ls->name(), *std::static_pointer_cast< QI::Lineshapes::NAME >(ls))); }
    if QI_SAVE( Gaussian )
    else if QI_SAVE( Splineshape )
    else { QI_FAIL("Unimplemented save for lineshape: " << ls->name()); }
    #undef QI_SAVE
}

void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::Lineshapes::Lineshape> &ls) {
    std::string type = ar.getNodeName();
    #define QI_LOAD( NAME ) \
        (type == #NAME ) { QI::Lineshapes::NAME l; ar(l); ls = std::make_shared< QI::Lineshapes::NAME >(l); }
    if QI_LOAD( Gaussian )
    else if QI_LOAD( Splineshape )
    else { QI_FAIL("Unimplemented load for lineshape: " << type); }
    #undef QI_LOAD
}

} // End namespace cereal

