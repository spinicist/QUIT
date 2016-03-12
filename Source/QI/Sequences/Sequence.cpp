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

#include "QI/Sequences/Sequence.h"
#include "QI/Util.h"

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
