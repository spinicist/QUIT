/*
 *  SpinEcho.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MultiEchoSequence.h"

namespace QI {

size_t MultiEchoSequence::size() const {
    return TE.rows();
}

Eigen::ArrayXcd MultiEchoSequence::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    return m->MultiEcho(p, TE, TR);
}

void MultiEchoSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(CEREAL_NVP(TR), CEREAL_NVP(TE1), CEREAL_NVP(ESP), CEREAL_NVP(ETL));
}

void MultiEchoSequence::load(cereal::JSONInputArchive &ar) {
    QI_SEQUENCE_LOAD( TR );
    QI_SEQUENCE_LOAD( TE1 );
    QI_SEQUENCE_LOAD( ESP );
    QI_SEQUENCE_LOAD( ETL );
    this->TE = Eigen::ArrayXd::LinSpaced(ETL, TE1, TE1+ESP*(ETL - 1));
}

} // End namespace QI
