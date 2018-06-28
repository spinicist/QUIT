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
#include "Macro.h"

namespace QI {
namespace Lineshapes {

Eigen::ArrayXd Gaussian::value(const Eigen::ArrayXd &df0, const double T2b) const {
    return sqrt(M_PI_2) * T2b * exp(-pow(2.0*M_PI*df0*T2b,2.0)/2.0);
}

void Gaussian::load(cereal::JSONInputArchive &ar) {
}

void Gaussian::save(cereal::JSONOutputArchive &ar) const {
}

} // End namespace Lineshapes
} // End namespace QI

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::Lineshapes::Lineshape> const &ls) {
    #define QI_SAVE( NAME ) \
        (ls->name() == #NAME ) { ar(cereal::make_nvp(ls->name(), *std::static_pointer_cast< QI::Lineshapes::NAME >(ls))); }
    if QI_SAVE( Gaussian )
    else { QI_FAIL("Unimplemented save for lineshape: " << ls->name()); }
    #undef QI_SAVE
}

void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::Lineshapes::Lineshape> &ls) {
    std::string type = ar.getNodeName();
    #define QI_LOAD( NAME ) \
        (type == #NAME ) { QI::Lineshapes::NAME l; ar(l); ls = std::make_shared< QI::Lineshapes::NAME >(l); }
    if QI_LOAD( Gaussian )
    else { QI_FAIL("Unimplemented load for lineshape: " << type); }
    #undef QI_LOAD
}

} // End namespace cereal

