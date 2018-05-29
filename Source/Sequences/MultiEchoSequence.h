/*
 *  MultiEchoSequence.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SPINECHO_H
#define SEQUENCES_SPINECHO_H

#include "SequenceBase.h"

namespace QI {

struct MultiEchoBase : SequenceBase {
    double TR;
    Eigen::ArrayXd TE;
    Eigen::Index size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<QI::Model> m, const Eigen::VectorXd &par) const override;
};
void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::MultiEchoBase> &sb);

struct MultiEchoSequence : MultiEchoBase {
    double TE1, ESP;
    int ETL;
    QI_SEQUENCE_DECLARE_NOSIG(MultiEcho);
};

struct MultiEchoFlexSequence : MultiEchoBase {
    QI_SEQUENCE_DECLARE_NOSIG(MultiEchoFlex);
};

} // End namespace QI

#endif // SEQUENCES_SPINECHO_H
