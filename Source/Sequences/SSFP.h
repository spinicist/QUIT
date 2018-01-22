/*
 *  SSFP.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SSFP_H
#define SEQUENCES_SSFP_H

#include "SequenceBase.h"
#include "Macro.h"

namespace QI {

struct SSFPBase : SequenceBase {
    double TR;
    Eigen::ArrayXd FA;

    size_t size() const override { return FA.rows(); }
};

struct SSFP : SSFPBase {
    Eigen::ArrayXd PhaseInc;

    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0) const override;
    template<typename Archive>
    void save(Archive &archive) const {
        archive(CEREAL_NVP(TR), CEREAL_NVP(FA), CEREAL_NVP(PhaseInc));
    }
    template<typename Archive>
    void load(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(FA), CEREAL_NVP(PhaseInc));
        if (FA.rows() != PhaseInc.rows()) {
            QI_FAIL("Number of flip-angles (FA) must equal number of phase-increments (Phase)");
        }
    }
};

struct SSFPEcho : SSFP {
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd signal_magnitude(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
};

struct SSFPFinite : SSFPBase {
    double Trf;
    Eigen::ArrayXd PhaseInc;

    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0) const override;
    template<typename Archive>
    void save(Archive &archive) const {
        archive(CEREAL_NVP(TR), CEREAL_NVP(Trf), CEREAL_NVP(FA), CEREAL_NVP(PhaseInc));
    }
    template<typename Archive>
    void load(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(Trf), CEREAL_NVP(FA), CEREAL_NVP(PhaseInc));
        if (FA.rows() != PhaseInc.rows()) {
            QI_FAIL("Number of flip-angles (FA) must equal number of phase-increments (Phase)");
        }
    }
};

struct SSFP_GS : SSFPBase {
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0) const override;
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(TR), CEREAL_NVP(FA));
    }
};

} // End namespace QI

CEREAL_REGISTER_TYPE(QI::SSFPBase);
CEREAL_REGISTER_POLYMORPHIC_RELATION(QI::SequenceBase, QI::SSFPBase);
CEREAL_REGISTER_TYPE(QI::SSFP);
CEREAL_REGISTER_POLYMORPHIC_RELATION(QI::SSFPBase, QI::SSFP);
CEREAL_REGISTER_TYPE(QI::SSFPEcho);
CEREAL_REGISTER_POLYMORPHIC_RELATION(QI::SSFP, QI::SSFPEcho);
CEREAL_REGISTER_TYPE(QI::SSFPFinite);
CEREAL_REGISTER_POLYMORPHIC_RELATION(QI::SSFPBase, QI::SSFPFinite);
CEREAL_REGISTER_TYPE(QI::SSFP_GS);
CEREAL_REGISTER_POLYMORPHIC_RELATION(QI::SSFPBase, QI::SSFP_GS);

#endif // SEQUENCES_SSFP_H
