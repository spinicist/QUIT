/*
 *  SpinEcho.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Sequences/SpinEcho.h"
#include "QI/Util.h"

using namespace std;
using namespace Eigen;

namespace QI {

/******************************************************************************
 * MultiEcho
 *****************************************************************************/
MultiEcho::MultiEcho() : SequenceBase() {}
MultiEcho::MultiEcho(const ArrayXd &te) :
	m_TE(te)
{}

MultiEcho::MultiEcho(std::istream& istr, const bool prompt) : SequenceBase() {
	double TE1;
	int NE;
    if (prompt) cout << "Enter TR: " << flush;
    QI::Read(istr, m_TR);
	if (prompt) cout << "Enter first echo-time: " << flush;
	QI::Read(istr, TE1);
	if (prompt) cout << "Enter echo spacing: " << flush;
	QI::Read(istr, m_ESP);
	if (prompt) cout << "Enter number of echos: " << flush;
	QI::Read(istr, NE);
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

} // End namespace QI
