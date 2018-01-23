/*
 *  MPRAGESequence.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_MPRAGE_H
#define SEQUENCES_MPRAGE_H

#include "SequenceBase.h"
#include "Macro.h"
#include "EigenCereal.h"

namespace QI {

struct MPRAGE : SequenceBase {
    double TR, FA, eta, TI, TD;
    int ETL, k0;
    std::string &name() const override { static std::string name = "MPRAGE"; return name; }
    size_t size() const override { return 1; }
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0 = 0.0) const override;
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(FA), CEREAL_NVP(eta), CEREAL_NVP(ETL), CEREAL_NVP(k0),
                CEREAL_NVP(TI), CEREAL_NVP(TD));
    }
};

struct MP2RAGE {
    double TR;
    int ETL;
    Eigen::ArrayXd FA;
    Eigen::ArrayXd TD;
    Eigen::ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
    template<typename Archive>
    void serialize(Archive &archive) {
        
        archive(CEREAL_NVP(TR), CEREAL_NVP(ETL), CEREAL_NVP(FA), CEREAL_NVP(TD));
    }
};

class MP3RAGE {
    double TR;
    int ETL;
    Eigen::ArrayXd FA;
    Eigen::ArrayXd TD;
    Eigen::ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(ETL), CEREAL_NVP(FA), CEREAL_NVP(TD));
    }
};

} // End namespace QI

#endif // SEQUENCES_MPRAGE_H
