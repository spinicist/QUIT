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

#define FA_PHASE_CHECK()\
    if (FA.rows() != PhaseInc.rows()) {\
        QI_FAIL("While reading " << this->name() << " number of phase increments " << PhaseInc.rows() <<\
                " did not match the number of flip angles " << FA.rows());\
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

SSFPSequence::SSFPSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
    PhaseInc = ArrayFromJSON(json["PhaseInc"], M_PI / 180);
    FA_PHASE_CHECK()
}

rapidjson::Value SSFPSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    json.AddMember("PhaseInc", ArrayToJSON(PhaseInc, a, 180 / M_PI), a);
    return json;
}

Eigen::ArrayXcd SSFPEchoSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFPEcho(p, FA, TR, PhaseInc);
}

Eigen::ArrayXd SSFPEchoSequence::signal_magnitude(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFPEchoMagnitude(p, FA, TR, PhaseInc);
}

SSFPEchoSequence::SSFPEchoSequence(const rapidjson::Value &json) :
    SSFPSequence(json)
{
}

rapidjson::Value SSFPEchoSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    json.AddMember("PhaseInc", ArrayToJSON(PhaseInc, a, 180 / M_PI), a);
    return json;
}

Eigen::ArrayXd SSFPFiniteSequence::weights(const double f0) const {
    Eigen::ArrayXd offset = PhaseInc + 2.*M_PI*f0*TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

Eigen::ArrayXcd SSFPFiniteSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFPFinite(p, FA, TR, Trf, PhaseInc);
}

SSFPFiniteSequence::SSFPFiniteSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    Trf = json["Trf"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
    PhaseInc = ArrayFromJSON(json["PhaseInc"], M_PI / 180);
    FA_PHASE_CHECK()
}

rapidjson::Value SSFPFiniteSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("Trf", Trf, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    json.AddMember("PhaseInc", ArrayToJSON(PhaseInc, a, 180 / M_PI), a);
    return json;
}

Eigen::ArrayXcd SSFPGSSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SSFP_GS(p, FA, TR);
}

SSFPGSSequence::SSFPGSSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
}

rapidjson::Value SSFPGSSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    return json;

}

Eigen::ArrayXcd SSFPEllipseSequence::signal(std::shared_ptr<Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not implemented");
}

Eigen::Index SSFPEllipseSequence::size() const {
    return FA.rows() * PhaseInc.rows();
}

SSFPEllipseSequence::SSFPEllipseSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
    PhaseInc = ArrayFromJSON(json["PhaseInc"], M_PI / 180);
}

rapidjson::Value SSFPEllipseSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    json.AddMember("PhaseInc", ArrayToJSON(PhaseInc, a, 180 / M_PI), a);
    return json;
}

Eigen::ArrayXcd SSFPMTSequence::signal(std::shared_ptr<Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not implemented");
}

Eigen::Index SSFPMTSequence::size() const {
    return FA.rows();
}

SSFPMTSequence::SSFPMTSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = ArrayFromJSON(json["TR"]);
    Trf = ArrayFromJSON(json["Trf"]);
    intB1 = ArrayFromJSON(json["intB1"]);
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
    PhaseInc = ArrayFromJSON(json["PhaseInc"], M_PI / 180);
    if ((TR.rows() != Trf.rows()) || (TR.rows() != intB1.rows()) || (TR.rows() != FA.rows())) {
        QI_FAIL("One on more parameters had differing lengths, TR had " << TR.rows() << ", Trf had " << Trf.rows() << ", intB1 had " << intB1.rows() << ", FA had " << FA.rows());
    }
}

rapidjson::Value SSFPMTSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", ArrayToJSON(TR, a), a);
    json.AddMember("Trf", ArrayToJSON(Trf, a), a);
    json.AddMember("intB1", ArrayToJSON(intB1, a), a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    json.AddMember("PhaseInc", ArrayToJSON(PhaseInc, a, 180 / M_PI), a);
    return json;
}

} // End namespace QI
