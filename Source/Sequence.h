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

#ifndef DESPOT_SEQUENCE
#define DESPOT_SEQUENCE

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Dense>

#include "SignalEquations.h"
#include "Model.h"

using namespace std;
using namespace Eigen;

class SequenceBase {
	protected:
		double m_TR = 0.;
		ArrayXd m_flip = ArrayXd::Zero(1);

	public:
		virtual ArrayXcd signal(const shared_ptr<Model> m, const VectorXd &p) const = 0;
		virtual size_t size() const = 0;
		virtual void write(ostream &os) const = 0;
		virtual string name() const = 0;
		virtual size_t count() const { return 1; }
		double TR() const { return m_TR; }
		void setTR(const double TR) { m_TR = TR; }
		const ArrayXd & flip() const { return m_flip; }
		void setFlip(const ArrayXd &f) { m_flip = f; }
};

ostream& operator<<(ostream& os, const SequenceBase& s);

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
};

class SteadyState : public SequenceBase {
	public:
		SteadyState();
		SteadyState(const ArrayXd &flip, const double TR);

		virtual size_t size() const override { return angles() * phases(); }
		virtual size_t angles() const { return m_flip.rows(); }
		virtual size_t phases() const { return 1; }
};

class SPGRSimple : public SteadyState {
	public:
		SPGRSimple(const ArrayXd &flip, const double TR);
		SPGRSimple(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR"; }
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
class MPRAGE : public SteadyState {
	public:
		ArrayXd m_TI;
		double m_TD;
		int m_N;
		MPRAGE() : SteadyState() {}
		MPRAGE(const ArrayXd &TI, const double TD, const double TR, const int N, const double flip);
		MPRAGE(const bool prompt);
		size_t size() const override { return m_TI.size(); }
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "MPRAGE"; }
};

// Special class for GE IRSPGR, for backwards compatibility
class IRSPGR : public MPRAGE {
	public:
		IRSPGR(const bool prompt);
		string name() const override { return "IRSPGR"; }
};

class SSFPSimple : public SteadyState {
	public:
		ArrayXd m_phases;
		SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases);
		SSFPSimple(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		size_t phases() const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP"; }

		bool isSymmetric() const;
		Array2d bandwidth() const;
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
		virtual size_t angles() const { return m_flip.rows(); }
		virtual size_t phases() const { return 1; }

		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "AFI"; }
};

class SequenceGroup : public SequenceBase {
private:
	vector<shared_ptr<SteadyState>> m_sequences;

public:
	SequenceGroup();
	void write(ostream &os) const override;
	string name() const override { return "Sequences"; }

	size_t count() const override;
	shared_ptr<SteadyState> sequence(const size_t i) const;
	vector<shared_ptr<SteadyState>> &sequences();

	size_t size() const override;
	ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
	//ArrayXcd loadSignals(vector<QUIT::MultiArray<complex<float>, 4>> &sigs, const size_t i, const size_t j, const size_t k, bool needsFlip = false) const;
	
	void addSequence(const shared_ptr<SteadyState> &seq);
};

#endif
