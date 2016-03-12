/*
 *  SteadyStateSequence.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_STEADYSTATE_H
#define SEQUENCES_STEADYSTATE_H

#include "QI/Sequences/SequenceBase.h"

class SteadyState : public SequenceBase {
    public:
        SteadyState();
        SteadyState(const ArrayXd &flip, const double TR);

        virtual size_t size() const override { return m_flip.rows(); }
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

class SSFP_GS : public SteadyState {
    public:
        SSFP_GS(const bool prompt);
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

#endif // SEQUENCES_STEADYSTATE_H