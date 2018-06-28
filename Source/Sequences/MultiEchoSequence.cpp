/*
 *  MultiEchoSequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MultiEchoSequence.h"
#include "CerealMacro.h"
#include "CerealEigen.h"

namespace QI {

/*
 * Base
 */
Eigen::Index MultiEchoBase::size() const {
    return TE.rows();
}

Eigen::ArrayXcd MultiEchoBase::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->MultiEcho(p, TE, TR);
}

void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::MultiEchoBase> &sb) {
    std::string seq_type = ar.getNodeName();
    #define QI_LOAD( NAME ) \
        (seq_type == #NAME ) { QI::NAME ## Sequence s; ar(s); sb = std::make_shared< QI::NAME ## Sequence >(s); }
    if QI_LOAD( MultiEcho )
    else if QI_LOAD( MultiEchoFlex )
    else { QI_FAIL("Not a Multi-Echo sequence type: " << seq_type); }
    #undef QI_LOAD
}

/*
 * Regularly spaced sequence
 */
void MultiEchoSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, TE1);
    QI_CSAVE(ar, ESP);
    QI_CSAVE(ar, ETL);
}

void MultiEchoSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, TE1);
    QI_CLOAD(ar, ESP);
    QI_CLOAD(ar, ETL);
    this->TE = Eigen::ArrayXd::LinSpaced(ETL, TE1, TE1+ESP*(ETL - 1));
}

/*
 * Irregularly spaced sequence
 */
void MultiEchoFlexSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, TE);
}

void MultiEchoFlexSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, TE);
}

} // End namespace QI
