/*
 *  Sequence.cpp
 *
 *  Created by Tobias Wood on 14/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Sequence.h"
#include "Util.h"

/******************************************************************************
 * MultiEcho
 *****************************************************************************/
MultiEcho::MultiEcho() : SequenceBase() {}
MultiEcho::MultiEcho(const ArrayXd &te) :
	m_TE(te)
{}

MultiEcho::MultiEcho(const bool prompt) : SequenceBase() {
	double TE1;
	int NE;
    if (prompt) cout << "Enter TR: " << flush;
    QI::Read(cin, m_TR);
	if (prompt) cout << "Enter first echo-time: " << flush;
	QI::Read(cin, TE1);
	if (prompt) cout << "Enter echo spacing: " << flush;
	QI::Read(cin, m_ESP);
	if (prompt) cout << "Enter number of echos: " << flush;
	QI::Read(cin, NE);
	m_TE.resize(NE);
	m_TE(0) = TE1;
	for (int i = 1; i < NE; i++) {
		m_TE(i) = m_TE(i-1) + m_ESP;
	}
}

ArrayXcd MultiEcho::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->MultiEcho(p, m_TE, m_TR);
}

void MultiEcho::write(ostream &os) const {
	os << "MultiEcho" << endl;
	os << "TE: " << m_TE.transpose() << endl;
}

SPGRSimple::SPGRSimple(const ArrayXd &flip, const double TR) :
	SteadyState(flip, TR)
{}
SPGRSimple::SPGRSimple(const bool prompt) : SteadyState() {
	if (prompt) cout << "Enter flip-angles (degrees): " << flush;
	QI::ReadArray(cin, m_flip);
	m_flip *= M_PI / 180.;
	if (prompt) cout << "Enter TR (seconds): " << flush;
	QI::Read(cin, m_TR);
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

SPGREcho::SPGREcho(const bool prompt) : SPGRSimple(prompt) {
    if (prompt) cout << "Enter TE (seconds): " << flush;
    QI::Read(cin, m_TE);
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
SPGRFinite::SPGRFinite(const bool prompt) : SPGRSimple(prompt) {
	if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
	QI::Read(cin, m_Trf);
	if (prompt) cout << "Enter TE (seconds): " << flush;
	QI::Read(cin, m_TE);
}

void SPGRFinite::write(ostream &os) const {
	os << "SPGR Finite" << endl;
	os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tTE: " << m_TE << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRFinite::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->SPGRFinite(p, m_flip, m_TR, m_Trf, m_TE);
}

SSFPSimple::SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phi) :
    SteadyState()
{
    m_TR = TR;
    m_flip = flip * M_PI / 180.;
    m_phi = phi * M_PI / 180.;
}

SSFPSimple::SSFPSimple(const bool prompt) :
	SteadyState()
{
	if (prompt) cout << "Enter flip-angles (degrees): " << flush;
    QI::ReadArray(cin, m_flip);
    if (prompt) cout << "Enter phase-increments (degrees): " << flush;
    QI::ReadArray(cin, m_phi);
    if (m_flip.size() != m_phi.size()) {
        QI_EXCEPTION("SSFP must have the same number of flip-angles and phase-increments.");
    }
	if (prompt) cout << "Enter TR (seconds): " << flush;
	QI::Read(cin, m_TR);
    m_flip *= M_PI / 180.;
    m_phi *= M_PI / 180.;
}

void SSFPSimple::write(ostream &os) const {
    os << name() << endl;
    os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
    os << "Phases: " << (m_phi * 180. / M_PI).transpose() << endl;
    os << "TR: " << m_TR << endl;
}

ArrayXd SSFPSimple::weights(const double f0) const {
    ArrayXd offset = m_phi + 2.*M_PI*f0*m_TR;
	ArrayXd weights = 0.75 * (offset / 2).sin().square();
	return weights;
}

ArrayXcd SSFPSimple::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFP(p, m_flip, m_TR, m_phi);
}
ArrayXcd SSFPEcho::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFPEcho(p, m_flip, m_TR, m_phi);
}
SSFPFinite::SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases) :
	SSFPSimple(flip, TR, phases), m_Trf(Trf)
{}
SSFPFinite::SSFPFinite(const bool prompt) :
	SSFPSimple(prompt)
{
	if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
	QI::Read(cin, m_Trf);
}

void SSFPFinite::write(ostream &os) const {
	os << "SSFP Finite" << endl;
    os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tPhases: " << (m_phi * 180. / M_PI).transpose() << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SSFPFinite::signal(shared_ptr<Model> m, const VectorXd &p) const {
    return m->SSFPFinite(p, m_flip, m_TR, m_Trf, m_phi);
}

SSFPEllipse::SSFPEllipse(const bool prompt) :
	SteadyState()
{
	if (prompt) cout << "Enter flip-angles (degrees): " << flush;
	QI::ReadArray(cin, m_flip);
	m_flip *= M_PI / 180.;
	if (prompt) cout << "Enter TR (seconds): " << flush;
	QI::Read(cin, m_TR);
}

void SSFPEllipse::write(ostream &os) const {
	os << "SSFP Ellipse" << endl;
	os << "TR: " << m_TR << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SSFPEllipse::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->SSFPEllipse(p, m_flip, m_TR);
}

AFI::AFI(const bool prompt) : SteadyState() {
	if (prompt) cout << "Enter flip-angle (degrees): " << flush;
	double inFlip;
	QI::Read(cin, inFlip);
	m_flip = ArrayXd::Ones(1) * inFlip * M_PI / 180.;
	if (prompt) cout << "Enter TR1 & TR2 (seconds): " << flush;
	ArrayXd temp;
	QI::ReadArray(cin, temp);
	if (temp.rows() != 2)
        QI_EXCEPTION("Must enter 2 TR values.");
	m_TR1 = temp[0]; m_TR2 = temp[1];
}

void AFI::write(ostream &os) const {
	os << name() << endl;
	os << "TR1: " << m_TR1 << " TR2: " << m_TR2 << endl;
	os << "Angle: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd AFI::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->AFI(p, m_flip[0], m_TR1, m_TR2);
}
