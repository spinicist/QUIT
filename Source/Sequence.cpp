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
 * SequenceBase
 *****************************************************************************/
ostream& operator<<(ostream& os, const SequenceBase& s) {
	s.write(os);
	return os;
}

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

/******************************************************************************
 * SteadyState
 *****************************************************************************/
SteadyState::SteadyState() : SequenceBase() {}
SteadyState::SteadyState(const ArrayXd &flip, const double TR) :
	SequenceBase()
{
	m_TR = TR;
	m_flip = flip;
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

MPRAGE::MPRAGE(const ArrayXd &TI, const ArrayXd &TD, const double TR, const int Nseg, const int Nk0, const double flip) :
    SteadyState(), m_TI(TI), m_TD(TD), m_Nseg(Nseg), m_Nk0(Nk0) {
	m_TR = TR;
	m_flip.resize(1); m_flip[0] = flip;
    if (m_TI.size() != m_TD.size()) {
        QI_EXCEPTION("Inversion and delay time arrays must match in size.");
    }
}

MPRAGE::MPRAGE(const bool prompt) : SteadyState() {
	if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
    double inFlip, TRseg;
	QI::Read(cin, inFlip);
	m_flip = ArrayXd::Ones(1) * inFlip * M_PI / 180.;
	if (prompt) cout << "Enter read-out TR (seconds): " << flush;
	QI::Read(cin, m_TR);
	if (prompt) cout << "Enter segment size: " << flush;
	QI::Read(cin, m_Nseg);
    if (prompt) cout << "Enter k0: " << flush;
    QI::Read(cin, m_Nk0);
	if (prompt) cout << "Enter inversion times (seconds): " << flush;
	QI::ReadArray(cin, m_TI);
    if (prompt) cout << "Enter relaxation delay times (seconds): " << flush;
    QI::ReadArray(cin, m_TD);
}

IRSPGR::IRSPGR(const bool prompt) : MPRAGE() {
	if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
	double inFlip;
	QI::Read(cin, inFlip);
	m_flip = ArrayXd::Ones(1) * inFlip * M_PI / 180.;
	if (prompt) cout << "Enter read-out TR (seconds): " << flush;
	QI::Read(cin, m_TR);

	int NPE2;
    if (prompt) cout << "Enter number of spatial locations (remember +4): ";
	QI::Read(cin, NPE2);
    if (NPE2 >= 64) {
        m_Nseg = NPE2 / 2;
    } else {
        m_Nseg = NPE2;
    }
    m_Nk0 = 0;
    
	if (prompt) cout << "Enter TIs (seconds): " << flush;
	QI::ReadArray(cin, m_TI);

    m_TD = ArrayXd::Zero(m_TI.size()); // For GE IR-SPGR the delay time is zero
}

ArrayXcd MPRAGE::signal(shared_ptr<Model> m, const VectorXd &par) const {
    return m->MPRAGE(par, m_flip[0], m_TR, m_Nseg, m_Nk0, m_TI, m_TD);
}

void MPRAGE::write(ostream &os) const {
	os << name() << endl;
    os << "TR: " << m_TR << "\tSegment Length: " << m_Nseg << "\tk-Zero: " << m_Nk0 << "\tAlpha: " << m_flip[0] * 180 / M_PI << endl;
	os << "TI: " << m_TI.transpose() << endl;
    os << "TD: " << m_TD.transpose() << endl;
}

ArrayXd MPRAGE::weights(const double f0) const {
    return ArrayXd::Ones(size()) * 1.0;
}

MP2RAGE::MP2RAGE(const Array3d &TD, const double TR, const int N, const Array2d flip) :
    SteadyState(), m_TD(TD), m_N(N) {
    m_TR = TR;
    m_flip = flip;
}

MP2RAGE::MP2RAGE(const bool prompt) : SteadyState() {
    if (prompt) cout << "Enter read-out flip-angles (degrees): " << flush;
    QI::ReadArray(cin, m_flip);
    if (m_flip.size() != 2) {
        QI_EXCEPTION("Must have 2 flip-angles for MP3RAGE");
    }
    m_flip *= M_PI / 180.;
    if (prompt) cout << "Enter read-out TR (seconds): " << flush;
    QI::Read(cin, m_TR);
    if (prompt) cout << "Enter segment size: " << flush;
    QI::Read(cin, m_N);
    if (prompt) cout << "Enter inversion times (seconds): " << flush;
    ArrayXd m_TI;
    QI::ReadArray(cin, m_TI);
    if (m_TI.size() != 2) {
        QI_EXCEPTION("Must specify 2 TI times for MP3RAGE");
    }
    // Assume centric for now
    m_TD[0] = m_TI[0];
    m_TD[1] = m_TI[1] - m_TD[0] - m_N*m_TR;
    if (prompt) cout << "Enter overall TR (seconds): " << flush;
    float TRseg;
    QI::Read(cin, TRseg);
    m_TD[2] = TRseg - m_TD[1] - m_N*m_TR;
}

ArrayXcd MP2RAGE::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP2RAGE(m_flip, m_TR, m_N, m_TD, M0, T1, B1, eta);
}

void MP2RAGE::write(ostream &os) const {
    os << name() << endl;
    os << "TR: " << m_TR << "\tN: " << m_N << "\tAlpha: " << m_flip[0] * 180 / M_PI << endl;
    os << "TD: " << m_TD.transpose() << "\tTS: " << (m_TD.sum() + 2*m_N*m_TR) << endl;
}

MP3RAGE::MP3RAGE(const Array4d &TD, const double TR, const int N, const Array3d flip) :
    SteadyState(), m_TD(TD), m_N(N) {
    m_TR = TR;
    m_flip = flip;
}

MP3RAGE::MP3RAGE(const bool prompt) : SteadyState() {
    if (prompt) cout << "Enter read-out flip-angles (degrees): " << flush;
    QI::ReadArray(cin, m_flip);
    if (m_flip.size() != 3) {
        QI_EXCEPTION("Must have 3 flip-angles for MP3RAGE");
    }
    m_flip *= M_PI / 180.;
    if (prompt) cout << "Enter read-out TR (seconds): " << flush;
    QI::Read(cin, m_TR);
    if (prompt) cout << "Enter segment size: " << flush;
    QI::Read(cin, m_N);
    if (prompt) cout << "Enter inversion times (seconds): " << flush;
    ArrayXd m_TI;
    QI::ReadArray(cin, m_TI);
    if (m_TI.size() != 3) {
        QI_EXCEPTION("Must specify 3 TI times for MP3RAGE");
    }
    // Assume centric for now
    m_TD[0] = m_TI[0];
    m_TD[1] = m_TI[1] - m_TD[0] - m_N*m_TR;
    m_TD[2] = m_TI[2] - m_TD[1] - m_N*m_TR;
    if (prompt) cout << "Enter overall TR (seconds): " << flush;
    float TRseg;
    QI::Read(cin, TRseg);
    m_TD[3] = TRseg - m_TD[2] - m_N*m_TR;
}

ArrayXcd MP3RAGE::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP3RAGE(m_flip, m_TR, m_N, m_TD, M0, T1, B1, eta);
}

void MP3RAGE::write(ostream &os) const {
    os << name() << endl;
    os << "TR: " << m_TR << "\tN: " << m_N << "\tAlpha: " << m_flip[0] * 180 / M_PI << endl;
    os << "TD: " << m_TD.transpose() << "\tTS: " << (m_TD.sum() + 3*m_N*m_TR) << endl;
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
    ArrayXd offset = m_phi + (M_PI * f0*m_TR);
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
/******************************************************************************
 * SequenceGroup Class
 *****************************************************************************/
SequenceGroup::SequenceGroup() : SequenceBase()
{}

void SequenceGroup::write(ostream &os) const {
	os << "Combined Sequence Count: " << m_sequences.size() << "\tCombined size: " << size() << endl;
	for (auto& s : m_sequences)
		os << *s;
}

size_t SequenceGroup::count() const {
	return m_sequences.size();
}

shared_ptr<SteadyState> SequenceGroup::sequence(const size_t i) const {
	return m_sequences.at(i);
}

vector<shared_ptr<SteadyState>> &SequenceGroup::sequences() {
	return m_sequences;
}

size_t SequenceGroup::size() const {
	size_t sz = 0;
	for (auto& sig : m_sequences)
		sz += sig->size();
	return sz;
}

ArrayXcd SequenceGroup::signal(shared_ptr<Model> m, const VectorXd &p) const {
	ArrayXcd result(size());
	size_t start = 0;
	for (auto &sig : m_sequences) {
		ArrayXcd thisResult = sig->signal(m, p);
		result.segment(start, sig->size()) = thisResult;
		start += sig->size();
	}
	return result;
}

ArrayXd SequenceGroup::weights(const double f0) const {
	ArrayXd weights(size());
	size_t start = 0;
	for (auto &sig : m_sequences) {
		weights.segment(start, sig->size()) = sig->weights(f0);
		start += sig->size();
	}
	return weights;
}

void SequenceGroup::addSequence(const shared_ptr<SteadyState> &seq) {
	m_sequences.push_back(seq);
}
