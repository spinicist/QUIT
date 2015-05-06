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
	if (prompt) cout << "Enter first echo-time: " << flush;
	QUITK::Read(cin, TE1);
	if (prompt) cout << "Enter echo spacing: " << flush;
	QUITK::Read(cin, m_ESP);
	if (prompt) cout << "Enter number of echos: " << flush;
	QUITK::Read(cin, NE);
	m_TE.resize(NE);
	m_TE(0) = TE1;
	for (int i = 1; i < NE; i++) {
		m_TE(i) = m_TE(i-1) + m_ESP;
	}
}

ArrayXcd MultiEcho::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->MultiEcho(p, m_TE);
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
	size_t nFlip;
	if (prompt) cout << "Enter number of flip-angles: " << flush;
	QUITK::Read(cin, nFlip);
	ArrayXd inAngles(nFlip);
	if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
	QUITK::ReadEigen(cin, inAngles);
	m_flip = inAngles * M_PI / 180.;
	if (prompt) cout << "Enter TR (seconds): " << flush;
	QUITK::Read(cin, m_TR);
}

void SPGRSimple::write(ostream &os) const {
	os << "SPGR Simple" << endl;
	os << "TR: " << m_TR << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRSimple::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->SPGR(p, m_flip, m_TR);
}

SPGRFinite::SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE) :
	SPGRSimple(flip, TR), m_Trf(Trf), m_TE(TE)
{}
SPGRFinite::SPGRFinite(const bool prompt) : SPGRSimple(prompt) {
	if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
	QUITK::Read(cin, m_Trf);
	if (prompt) cout << "Enter TE (seconds): " << flush;
	QUITK::Read(cin, m_TE);
}

void SPGRFinite::write(ostream &os) const {
	os << "SPGR Finite" << endl;
	os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tTE: " << m_TE << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRFinite::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->SPGRFinite(p, m_flip, m_TR, m_Trf, m_TE);
}

MPRAGE::MPRAGE(const ArrayXd &TI, const double TD, const double TR, const int N, const double flip) :
	SteadyState(), m_TI(TI), m_TD(TD), m_N(N) {
	m_TR = TR;
	m_flip.resize(1); m_flip[0] = flip;
}

MPRAGE::MPRAGE(const bool prompt) : SteadyState() {
	if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
	ArrayXd inFlip(1);
	QUITK::ReadEigen(cin, inFlip);
	m_flip = inFlip * M_PI / 180.;
	if (prompt) cout << "Enter read-out TR (seconds): " << flush;
	QUITK::Read(cin, m_TR);
	if (prompt) cout << "Enter segment size: " << flush;
	QUITK::Read(cin, m_N);
	size_t nTI;
	if (prompt) cout << "Enter number of inversion times: " << flush;
	QUITK::Read(cin, nTI);
	m_TI.resize(nTI);
	if (prompt) cout << "Enter " << m_TI.size() << " inversion times (seconds): " << flush;
	QUITK::ReadEigen(cin, m_TI);
	if (prompt) cout << "Enter delay time (seconds): " << flush;
	QUITK::Read(cin, m_TD);
}

IRSPGR::IRSPGR(const bool prompt) : MPRAGE() {
	if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
	ArrayXd inFlip(1);
	QUITK::ReadEigen(cin, inFlip);
	m_flip = inFlip * M_PI / 180.;
	if (prompt) cout << "Enter read-out TR (seconds): " << flush;
	QUITK::Read(cin, m_TR);

	int NPE2;
	if (prompt) cout << "Enter original number of slices (PE2):";
	QUITK::Read(cin, NPE2);
	m_N = (NPE2 / 2) + 2;

	int nTI;
	if (prompt) cout << "Enter number of TIs: " << flush;
	QUITK::Read(cin, nTI);
	m_TI.resize(nTI);
	if (prompt) cout << "Enter " << m_TI.size() << " TIs (seconds): " << flush;
	QUITK::ReadEigen(cin, m_TI);
}

ArrayXcd MPRAGE::signal(shared_ptr<Model> m, const VectorXd &par) const {
	return m->MPRAGE(par, m_flip[0], m_TR, m_N, m_TI, m_TD);
}

void MPRAGE::write(ostream &os) const {
	os << name() << endl;
	os << "TR: " << m_TR << "\tN: " << m_N << "\tAlpha: " << m_flip[0] * 180 / M_PI << "\tTD: " << m_TD << endl;
	os << "TI: " << m_TI.transpose() << endl;
	os << "TS: " << (m_TI + m_N*m_TR + m_TD).transpose() << endl;
}

SSFPSimple::SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases) :
	SteadyState(flip, TR), m_phases(phases)
{}
SSFPSimple::SSFPSimple(const bool prompt) :
	SteadyState()
{
	size_t nFlip;
	if (prompt) cout << "Enter number of flip-angles: " << flush;
	QUITK::Read(cin, nFlip);
	ArrayXd inAngles(nFlip);
	if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
	QUITK::ReadEigen(cin, inAngles);
	m_flip = inAngles * M_PI / 180.;
	size_t nPhases;
	if (prompt) cout << "Enter number of phase-cycles: " << flush;
	QUITK::Read(cin, nPhases);
	ArrayXd inPhases(nPhases);
	if (prompt) cout << "Enter " << inPhases.size() << " phase-cycles (degrees): " << flush;
	QUITK::ReadEigen(cin, inPhases);
	m_phases = inPhases * M_PI / 180.;
	if (prompt) cout << "Enter TR (seconds): " << flush;
	QUITK::Read(cin, m_TR);
}

void SSFPSimple::write(ostream &os) const {
	os << "SSFP Simple" << endl;
	os << "TR: " << m_TR << "\tPhases: " << (m_phases * 180. / M_PI).transpose() << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

size_t SSFPSimple::phases() const { return m_phases.rows(); }
bool SSFPSimple::isSymmetric() const {
	bool sym = true;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		if (!((abs(m_phases(i) - M_PI) <= (M_PI * numeric_limits<double>::epsilon())) ||
		      (abs(m_phases(i) - 0.) <= numeric_limits<double>::epsilon()))) {
			sym = false;
			break; // Don't need to bother checking other values
		}
	}
	return sym;
}

Array2d SSFPSimple::bandwidth() const {
	// For now, assume that if people have 4 phase-cycles they are sensible enough
	// to space them over 2*pi instead of pi.

	Array2d bw = Array2d::Zero();
	bw(1) = 0.5/m_TR;
	if (m_phases.rows() <= 2) {
		if (!isSymmetric())
			bw -= bw(1) / 2;
	} else {
		if (!isSymmetric())
			bw(0) = -bw(1);
	}
	return bw;
}

ArrayXcd SSFPSimple::signal(shared_ptr<Model> m, const VectorXd &p) const {
	ArrayXcd s(size());
	ArrayXcd::Index start = 0;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		s.segment(start, m_flip.rows()) = m->SSFP(p, m_flip, m_TR, m_phases(i));
		start += m_flip.rows();
	}
	return s;
}

SSFPFinite::SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases) :
	SSFPSimple(flip, TR, phases), m_Trf(Trf)
{}
SSFPFinite::SSFPFinite(const bool prompt) :
	SSFPSimple(prompt)
{
	if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
	QUITK::Read(cin, m_Trf);
}

void SSFPFinite::write(ostream &os) const {
	os << "SSFP Finite" << endl;
	os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tPhases: " << (m_phases * 180. / M_PI).transpose() << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SSFPFinite::signal(shared_ptr<Model> m, const VectorXd &p) const {
	ArrayXcd s(size());
	ArrayXcd::Index start = 0;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		s.segment(start, m_flip.rows()) = m->SSFPFinite(p, m_flip, m_TR, m_Trf, m_phases(i));
		start += m_flip.rows();
	}
	return s;
}

SSFPEllipse::SSFPEllipse(const bool prompt) :
	SteadyState()
{
	size_t nFlip;
	if (prompt) cout << "Enter number of flip-angles: " << flush;
	QUITK::Read(cin, nFlip);
	ArrayXd inAngles(nFlip);
	if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
	QUITK::ReadEigen(cin, inAngles);
	m_flip = inAngles * M_PI / 180.;
	if (prompt) cout << "Enter TR (seconds): " << flush;
	QUITK::Read(cin, m_TR);
}

void SSFPEllipse::write(ostream &os) const {
	os << "SSFP Ellipse" << endl;
	os << "TR: " << m_TR << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SSFPEllipse::signal(shared_ptr<Model> m, const VectorXd &p) const {
	return m->SSFPEllipse(p, m_flip, m_TR);
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

/*
ArrayXcd SequenceGroup::loadSignals(vector<QUITK::MultiArray<complex<float>, 4>> &sigs,
                                const size_t i, const size_t j, const size_t k,
                                const bool flip) const {
	ArrayXcd signal(size());
	size_t start = 0;
	for (size_t s = 0; s < m_sequences.size(); s++) {
		ArrayXcd thisSig = sigs.at(s).slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<complex<double>>();
		if (flip) {
			ArrayXXcd flipped = Map<ArrayXXcd>(thisSig.data(), m_sequences.at(s)->phases(), m_sequences.at(s)->angles()).transpose();
			thisSig = Map<ArrayXcd>(flipped.data(), thisSig.rows(), 1);
		}
		signal.segment(start, thisSig.rows()) = thisSig;
		start += thisSig.rows();
	}
	return signal;
}*/

void SequenceGroup::addSequence(const shared_ptr<SteadyState> &seq) {
	m_sequences.push_back(seq);
}
