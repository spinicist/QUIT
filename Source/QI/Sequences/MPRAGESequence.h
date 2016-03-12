/*
 *  MPRAGESequence.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_MPRAGE_H
#define SEQUENCES_MPRAGE_H

#include "QI/Sequences/SequenceBase.h"

namespace QI {

class MPRAGE : public SequenceBase {
    public:
        ArrayXd m_TI, m_TD;
        int m_Nseg, m_Nk0;
        MPRAGE() : SequenceBase() {}
        MPRAGE(const ArrayXd &TI, const ArrayXd &TD, const double TR, const int Nseg, const int Nk0, const double flip);
        MPRAGE(const bool prompt);
        size_t size() const override { return m_TI.size(); }
        ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
        void write(ostream &os) const override;
        string name() const override { return "MPRAGE"; }
        ArrayXd weights(const double f0 = 0.0) const override;
};
// Special class for GE IRSPGR, for backwards compatibility
class IRSPGR : public MPRAGE {
    public:
        IRSPGR(const bool prompt);
        string name() const override { return "IRSPGR"; }
};

class MP2RAGE : public SequenceBase {
    public:
        Array3d m_TD;
        int m_N;
        MP2RAGE() : SequenceBase() {}
        MP2RAGE(const Array3d &TD, const double TR, const int N, const Array2d flip);
        MP2RAGE(const bool prompt);
        size_t size() const override { return 3; }
        ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override { QI_EXCEPTION("Not implemented"); }
        ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
        void write(ostream &os) const override;
        string name() const override { return "MP3RAGE"; }
};

class MP3RAGE : public SequenceBase {
    public:
        Array4d m_TD;
        int m_N;
        MP3RAGE() : SequenceBase() {}
        MP3RAGE(const Array4d &TD, const double TR, const int N, const Array3d flip);
        MP3RAGE(const bool prompt);
        size_t size() const override { return 3; }
        ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override { QI_EXCEPTION("Not implemented"); }
        ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
        void write(ostream &os) const override;
        string name() const override { return "MP3RAGE"; }
};

} // End namespace QI

#endif // SEQUENCES_MPRAGE_H
