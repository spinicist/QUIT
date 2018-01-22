/*
 *  SpinEcho.h
 *
 *  Copyright (c) 2016 Tobias Wood.
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

struct MultiEcho : SequenceBase {
    double TR, ESP;
    int NEcho;
    Eigen::ArrayXd TE;

    size_t size() const override { return NEcho; }
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    template<typename Archive>
    void save(Archive &archive) const {
        archive(CEREAL_NVP(TR), CEREAL_NVP(ESP), CEREAL_NVP(NEcho));
    }
    template<typename Archive>
    void load(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(ESP), CEREAL_NVP(NEcho));
        this->TE = Eigen::ArrayXd::LinSpaced(NEcho, ESP, ESP*NEcho);
    }

};

} // End namespace QI

CEREAL_REGISTER_TYPE(QI::MultiEcho);
CEREAL_REGISTER_POLYMORPHIC_RELATION(QI::SequenceBase, QI::MultiEcho);

#endif // SEQUENCES_SPINECHO_H
