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
        ArrayXd m_flip;

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

        ArrayXd TE() const { return m_TE; }
        void setTE(const ArrayXd &TE) { m_TE = TE; }
};

class SteadyState : public SequenceBase {
	public:
		SteadyState();
		SteadyState(const ArrayXd &flip, const double TR);

        virtual size_t size() const override { return m_flip.rows(); }
		virtual ArrayXd weights(double f0) const { return ArrayXd::Ones(size()); }
};

class SPGRSimple : public SteadyState {
	public:
		SPGRSimple(const ArrayXd &flip, const double TR);
		SPGRSimple(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR"; }
		ArrayXd weights(const double f0) const override;
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
class MPRAGE : public SteadyState {
	public:
        ArrayXd m_TI, m_TRseg;
		int m_N;
		MPRAGE() : SteadyState() {}
        MPRAGE(const ArrayXd &TI, const double TRseg, const double TR, const int N, const double flip);
		MPRAGE(const bool prompt);
		size_t size() const override { return m_TI.size(); }
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "MPRAGE"; }
};
class MP2RAGE : public SteadyState {
    public:
        Array3d m_TD;
        int m_N;
        MP2RAGE() : SteadyState() {}
        MP2RAGE(const Array3d &TD, const double TR, const int N, const Array2d flip);
        MP2RAGE(const bool prompt);
        size_t size() const override { return 3; }
        ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override { throw(runtime_error(string(__PRETTY_FUNCTION__) + " unimplemented")); }
        ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
        void write(ostream &os) const override;
        string name() const override { return "MP3RAGE"; }
};
class MP3RAGE : public SteadyState {
    public:
        Array4d m_TD;
        int m_N;
        MP3RAGE() : SteadyState() {}
        MP3RAGE(const Array4d &TD, const double TR, const int N, const Array3d flip);
        MP3RAGE(const bool prompt);
        size_t size() const override { return 3; }
        ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override { throw(runtime_error(string(__PRETTY_FUNCTION__) + " unimplemented")); }
        ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
        void write(ostream &os) const override;
        string name() const override { return "MP3RAGE"; }
};
// Special class for GE IRSPGR, for backwards compatibility
class IRSPGR : public MPRAGE {
	public:
		IRSPGR(const bool prompt);
		string name() const override { return "IRSPGR"; }
};

class SSFPSimple : public SteadyState {
    private:
        size_t m_nphi;

	public:
        ArrayXd m_phi;
		SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases);
		SSFPSimple(const bool prompt);
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP"; }

		bool isSymmetric() const;
        virtual double bwMult() const;
        size_t phases() { return m_nphi; }
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

        virtual double bwMult() const override;
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
	ArrayXd weights(const double f0) const;
	
	void addSequence(const shared_ptr<SteadyState> &seq);
};

#endif
