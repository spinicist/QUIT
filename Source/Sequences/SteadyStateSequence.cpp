/*
 *  SteadyStateSequence.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SteadyStateSequence.h"
#include "IO.h"

using namespace std;
using namespace Eigen;

namespace QI {

/******************************************************************************
 * SteadyState
 *****************************************************************************/
 
SteadyState::SteadyState() : SequenceBase() {}

SteadyState::SteadyState(const ArrayXd &flip, const double TR) :
    SequenceBase(flip, TR)
{}

/******************************************************************************
 * SPGR
 *****************************************************************************/

SPGRSimple::SPGRSimple(const ArrayXd &flip, const double TR) :
    SteadyState(flip, TR)
{}
SPGRSimple::SPGRSimple(std::istream &istr, const bool prompt) : SteadyState() {
    if (prompt) cout << "Enter flip-angles (degrees): " << flush;
    QI::ReadArray(istr, m_flip);
    m_flip *= M_PI / 180.;
    if (prompt) cout << "Enter TR (seconds): " << flush;
    QI::Read(istr, m_TR);
}

void SPGRSimple::write(ostream &os) const {
    os << "SPGR Simple" << endl;
    os << "TR: " << m_TR << endl;
    os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRSimple::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SPGR(p, m_flip, m_TR);
}

ArrayXd SPGRSimple::weights(const double f0) const {
    // Weight SPGR images higher than SSFP
    return ArrayXd::Ones(size());
}

SPGREcho::SPGREcho(const ArrayXd &flip, const double TR, const double TE) :
    SPGRSimple(flip, TR), m_TE(TE)
{}

SPGREcho::SPGREcho(std::istream &istr, const bool prompt) : SPGRSimple(istr, prompt) {
    if (prompt) cout << "Enter TE (seconds): " << flush;
    QI::Read(istr, m_TE);
}

void SPGREcho::write(ostream &os) const {
    os << "SPGR Echo" << endl;
    os << "TR: " << m_TR << "\tTE: " << m_TE << endl;
    os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGREcho::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SPGREcho(p, m_flip, m_TR, m_TE);
}

SPGRFinite::SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE) :
    SPGRSimple(flip, TR), m_Trf(Trf), m_TE(TE)
{}
SPGRFinite::SPGRFinite(std::istream &istr, const bool prompt) : SPGRSimple(istr, prompt) {
    if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
    QI::Read(istr, m_Trf);
    if (prompt) cout << "Enter TE (seconds): " << flush;
    QI::Read(istr, m_TE);
}

void SPGRFinite::write(ostream &os) const {
    os << "SPGR Finite" << endl;
    os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tTE: " << m_TE << endl;
    os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRFinite::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SPGRFinite(p, m_flip, m_TR, m_Trf, m_TE);
}

/******************************************************************************
 * SSFP
 *****************************************************************************/
SSFPSimple::SSFPSimple() : SteadyState() {}

SSFPSimple::SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phi) :
    SteadyState(flip, TR),
    m_phi(phi * M_PI / 180.)
{
    setupAll();
}

SSFPSimple::SSFPSimple(std::istream &istr, const bool prompt) : SteadyState() {
    if (prompt) cout << "Enter flip-angles (degrees): " << flush;
    QI::ReadArray(istr, m_flip);
    if (prompt) cout << "Enter phase-increments (degrees): " << flush;
    QI::ReadArray(istr, m_phi);
    if (prompt) cout << "Enter TR (seconds): " << flush;
    QI::Read(istr, m_TR);
    m_flip *= M_PI / 180.;
    m_phi *= M_PI / 180.;
    setupAll();
}

void SSFPSimple::setupAll() {
    m_allFlip = m_flip.replicate(m_phi.rows(), 1);
    m_allPhi = ArrayXd::Zero(m_allFlip.size());
    int start = 0;
    for (int i = 0; i < m_phi.size(); i++) {
        m_allPhi.segment(start, m_flip.size()).setConstant(m_phi[i]);
        start += m_flip.size();
    }
}

void SSFPSimple::write(ostream &os) const {
    os << name() << endl;
    os << "Angles: " << (flip() * 180. / M_PI).transpose() << endl;
    os << "Phases: " << (phase_incs() * 180. / M_PI).transpose() << endl;
    os << "TR:     " << m_TR << endl;
}

ArrayXd SSFPSimple::weights(const double f0) const {
    ArrayXd offset = m_allPhi + 2.*M_PI*f0*m_TR;
    ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

ArrayXcd SSFPSimple::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFP(p, m_allFlip, m_TR, m_allPhi);
}

SSFPEcho::SSFPEcho() : SSFPSimple() {}
SSFPEcho::SSFPEcho(std::istream &istr, const bool prompt) : SSFPSimple(istr, prompt) {}

ArrayXcd SSFPEcho::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFPEcho(p, m_allFlip, m_TR, m_allPhi);
}

ArrayXd  SSFPEcho::signal_magnitude(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFPEchoMagnitude(p, m_allFlip, m_TR, m_allPhi);
}

SSFPEchoFlex::SSFPEchoFlex(std::istream &istr, const bool prompt) : SSFPEcho() {
    if (prompt) cout << "Enter flip-angles (degrees): " << flush;
    QI::ReadArray(istr, m_allFlip);
    if (prompt) cout << "Enter phase-increments (degrees): " << flush;
    QI::ReadArray(istr, m_allPhi);
    if (m_allFlip.size() != m_allPhi.size()) {
        QI_EXCEPTION("SSFP must have the same number of flip-angles and phase-increments." << std::endl <<
                     "Input flip-angles: " << m_allFlip.transpose() << std::endl <<
                     "Input phase-incs:  " << m_allPhi.transpose());
    }
    if (prompt) cout << "Enter TR (seconds): " << flush;
    QI::Read(istr, m_TR);
    m_allFlip *= M_PI / 180.;
    m_allPhi *= M_PI / 180.;
    m_flip = m_allFlip;
    m_phi = m_allPhi;
}

SSFPFinite::SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases) :
    SSFPSimple(flip, TR, phases),
    m_Trf(Trf)
{}

SSFPFinite::SSFPFinite(std::istream &istr, const bool prompt) : SSFPSimple(istr, prompt)
{
    if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
    QI::Read(istr, m_Trf);
}

void SSFPFinite::write(ostream &os) const {
    os << name() << endl;
    os << "Angles: " << (flip() * 180. / M_PI).transpose() << endl;
    os << "Phases: " << (phase_incs() * 180. / M_PI).transpose() << endl;
    os << "TR:     " << m_TR << "\tTrf: " << m_Trf << endl;
}

ArrayXcd SSFPFinite::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFPFinite(p, m_allFlip, m_TR, m_Trf, m_allPhi);
}

SSFP_GS::SSFP_GS(std::istream &istr, const bool prompt) : SteadyState() {
    if (prompt) cout << "Enter flip-angles (degrees): " << flush;
    QI::ReadArray(istr, m_flip);
    m_flip *= M_PI / 180.;
    if (prompt) cout << "Enter TR (seconds): " << flush;
    QI::Read(istr, m_TR);
}

void SSFP_GS::write(ostream &os) const {
    os << "SSFP Geometric Solution" << endl;
    os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
    os << "TR:     " << m_TR << endl;
}

ArrayXcd SSFP_GS::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFP_GS(p, m_flip, m_TR);
}

/******************************************************************************
 * AFI
 *****************************************************************************/

AFI::AFI(std::istream &istr, const bool prompt) : SteadyState() {
    if (prompt) cout << "Enter flip-angle (degrees): " << flush;
    double inFlip;
    QI::Read(istr, inFlip);
    m_flip = ArrayXd::Ones(1) * inFlip * M_PI / 180.;
    if (prompt) cout << "Enter TR1 & TR2 (seconds): " << flush;
    ArrayXd temp;
    QI::ReadArray(istr, temp);
    if (temp.rows() != 2)
        QI_EXCEPTION("Must enter 2 TR values.");
    m_TR1 = temp[0]; m_TR2 = temp[1];
}

void AFI::write(ostream &os) const {
    os << name() << endl;
    os << "Angle: " << (m_flip * 180. / M_PI).transpose() << endl;
    os << "TR1: " << m_TR1 << " TR2: " << m_TR2 << endl;
}

ArrayXcd AFI::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->AFI(p, m_flip[0], m_TR1, m_TR2);
}

} // End namespace QI
