/*
 *  Sequence.h
 *
 *  Created by Tobias Wood on 14/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SEQUENCE_H
#define SEQUENCES_SEQUENCE_H

#include "Sequences/SequenceBase.h"
#include "Sequences/SteadyStateSequence.h"

class MultiEcho : public SequenceBase {
	public:
		double m_ESP;
		ArrayXd m_TE;

		MultiEcho();
		MultiEcho(const ArrayXd &te);
		MultiEcho(const bool prompt);

		size_t size() const override { return m_TE.rows(); }
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "MultiEcho"; }

        ArrayXd TE() const { return m_TE; }
        void setTE(const ArrayXd &TE) { m_TE = TE; }
};

class SPGRSimple : public SteadyState {
	public:
		SPGRSimple(const ArrayXd &flip, const double TR);
		SPGRSimple(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR"; }
		ArrayXd weights(const double f0 = 0.0) const override;
};
class SPGREcho : public SPGRSimple {
public:
    double m_TE;
    SPGREcho(const ArrayXd &flip, const double TR, const double TE);
    SPGREcho(const bool prompt);
    ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
    void write(ostream &os) const override;
    string name() const override { return "SPGR_Echo"; }
};
class SPGRFinite : public SPGRSimple {
	public:
		double m_Trf, m_TE;
		SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE);
		SPGRFinite(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR_Finite"; }
};

class SSFPSimple : public SteadyState {
    protected:
        ArrayXd m_phi;
    
    public:
		SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases);
		SSFPSimple(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP"; }

		ArrayXd weights(const double f0) const override;
};
class SSFPEcho : public SSFPSimple {
public:
    SSFPEcho(const bool prompt) : SSFPSimple(prompt) {}

    ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
    string name() const override { return "SSFPEcho"; }
};
class SSFPFinite : public SSFPSimple {
	public:
		double m_Trf;
		SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases);
		SSFPFinite(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Finite"; }
};

class SSFPEllipse : public SteadyState {
	public:
		SSFPEllipse(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Ellipse"; }
};

class AFI : public SteadyState {
protected:
	double m_TR1, m_TR2;
	public:
		AFI(const bool prompt);

		virtual size_t size() const override { return 2; }
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "AFI"; }
};

#endif // SEQUENCES_SEQUENCE_H
