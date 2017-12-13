/*
 *  SpinEcho.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SEQUENCE_H
#define SEQUENCES_SEQUENCE_H

#include "SequenceBase.h"

namespace QI {

class MultiEcho : public SequenceBase {
	public:
		double m_ESP;
        Eigen::ArrayXd m_TE;

		MultiEcho();
        MultiEcho(const Eigen::ArrayXd &te);
		MultiEcho(std::istream &istr, const bool prompt);

		size_t size() const override { return m_TE.rows(); }
        Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
        void write(std::ostream &os) const override;
        std::string name() const override { return "MultiEcho"; }

        Eigen::ArrayXd TE() const { return m_TE; }
        void setTE(const Eigen::ArrayXd &TE) { m_TE = TE; }
};

} // End namespace QI

#endif // SEQUENCES_SEQUENCE_H
