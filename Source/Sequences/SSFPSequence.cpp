/*
 *  SSFP.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SSFPSequence.h"
#include "CerealMacro.h"

#define FA_PHASE_CHECK()\
    if (FA.rows() != PhaseInc.rows()) {\
        QI_FAIL("While reading " << this->name() << " number of phase increments " << PhaseInc.rows() <<\
                " did not match the number of flip angles " << FA.rows() );\
    }

namespace QI {

Eigen::Index SSFPBase::size() const {
    return FA.rows();
}

Eigen::ArrayXd SSFPSequence::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFPSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFP(p, FA, TR, PhaseInc);
}

void SSFPSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD_DEGREES(ar, FA);
    QI_CLOAD_DEGREES(ar, PhaseInc);
    FA_PHASE_CHECK()
}

void SSFPSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE_DEGREES(ar, FA);
    QI_CSAVE_DEGREES(ar, PhaseInc);
}

Eigen::ArrayXcd SSFPEchoSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFPEcho(p, FA, TR, PhaseInc);
}

Eigen::ArrayXd SSFPEchoSequence::signal_magnitude(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFPEchoMagnitude(p, FA, TR, PhaseInc);
}

void SSFPEchoSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD_DEGREES(ar, FA);
    QI_CLOAD_DEGREES(ar, PhaseInc);
    FA_PHASE_CHECK()
}

void SSFPEchoSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE_DEGREES(ar, FA);
    QI_CSAVE_DEGREES(ar, PhaseInc);
}

Eigen::ArrayXd SSFPFiniteSequence::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFPFiniteSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFPFinite(p, FA, TR, Trf, PhaseInc);
}

void SSFPFiniteSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, Trf);
    QI_CLOAD_DEGREES(ar, FA);
    QI_CLOAD_DEGREES(ar, PhaseInc);
    FA_PHASE_CHECK()
}

void SSFPFiniteSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, Trf);
    QI_CSAVE_DEGREES(ar, FA);
    QI_CSAVE_DEGREES(ar, PhaseInc);
}

Eigen::ArrayXcd SSFPGSSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFP_GS(p, FA, TR);
}

void SSFPGSSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD_DEGREES(ar, FA);
}

void SSFPGSSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE_DEGREES(ar, FA);
}

Eigen::ArrayXcd SSFPEllipseSequence::signal(std::shared_ptr<Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not implemented");
}

Eigen::Index SSFPEllipseSequence::size() const {
    return FA.rows() * PhaseInc.rows();
}

void SSFPEllipseSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD_DEGREES(ar, FA);
    QI_CLOAD_DEGREES(ar, PhaseInc);
}

void SSFPEllipseSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE_DEGREES(ar, FA);
    QI_CSAVE_DEGREES(ar, PhaseInc);
}

Eigen::ArrayXcd SSFPMTSequence::signal(std::shared_ptr<Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not implemented");
}

Eigen::Index SSFPMTSequence::size() const {
    return FA.rows();
}

void SSFPMTSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, Trf);
    QI_CLOAD(ar, intB1);
    QI_CLOAD_DEGREES(ar, FA);
    QI_CLOAD_DEGREES(ar, PhaseInc);
    if ((TR.rows() != Trf.rows()) || (TR.rows() != intB1.rows()) || (TR.rows() != FA.rows())) {
        QI_FAIL("One on more parameters had differing lengths, TR had " << TR.rows() << ", Trf had " << Trf.rows() << ", intB1 had " << intB1.rows() << ", FA had " << FA.rows());
    }
}

void SSFPMTSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, Trf);
    QI_CSAVE(ar, intB1);
    QI_CSAVE_DEGREES(ar, FA);
    QI_CSAVE_DEGREES(ar, PhaseInc);
}

} // End namespace QI
