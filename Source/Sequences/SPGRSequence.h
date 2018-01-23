/*
 *  SPGR.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR / FLASH / FFE Sequences
 *
 */

#ifndef SEQUENCES_SPGR_H
#define SEQUENCES_SPGR_H

#include "SequenceBase.h"
#include "cereal/cereal.hpp"
#include "EigenCereal.h"

namespace QI {

struct SPGRSequence : SequenceBase {
    double TR;
    Eigen::ArrayXd FA;

    std::string &name() const override { static std::string name = "SPGR"; return name; }
    size_t size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(FA));
    }
};

struct SPGREchoSequence : SequenceBase {
    double TR, TE;
    Eigen::ArrayXd FA;

    std::string &name() const override { static std::string name = "SPGREcho"; return name; }
    size_t size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(TE), CEREAL_NVP(FA));
    }
};

struct SPGRFiniteSequence : SequenceBase {
    double TR, TE, Trf;
    Eigen::ArrayXd FA;

    std::string &name() const override { static std::string name = "SPGRFinite"; return name; }
    size_t size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(TE), CEREAL_NVP(Trf), CEREAL_NVP(FA));
    }
};

} // End namespace QI

#endif // SEQUENCES_SPGR_H
