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
        Eigen::ArrayXd m_TI, m_TD;
        double m_eta;
        int m_Nseg, m_Nk0;
        MPRAGE() : SequenceBase() {}
        MPRAGE(const Eigen::ArrayXd &TI, const Eigen::ArrayXd &TD, const double TR, const int Nseg, const int Nk0, const double flip, const double eta);
        MPRAGE(std::istream &istr, const bool prompt);
        size_t size() const override { return m_TI.size(); }
        Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
        void write(std::ostream &os) const override;
        std::string name() const override { return "MPRAGE"; }
        Eigen::ArrayXd weights(const double f0 = 0.0) const override;
};
// Special class for GE IRSPGR, for backwards compatibility
class IRSPGR : public MPRAGE {
    public:
        IRSPGR(std::istream &istr, const bool prompt);
        std::string name() const override { return "IRSPGR"; }
};

class MP2RAGE : public SequenceBase {
    public:
        Eigen::Array3d m_TD;
        int m_N;
        MP2RAGE() : SequenceBase() {}
        MP2RAGE(const Eigen::Array3d &TD, const double TR, const int N, const Eigen::Array2d flip);
        MP2RAGE(std::istream &istr, const bool prompt);
        size_t size() const override { return 3; }
        Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override { QI_EXCEPTION("Not implemented"); }
        Eigen::ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
        void write(std::ostream &os) const override;
        std::string name() const override { return "MP3RAGE"; }
};

class MP3RAGE : public SequenceBase {
    public:
        Eigen::Array4d m_TD;
        int m_N;
        MP3RAGE() : SequenceBase() {}
        MP3RAGE(const Eigen::Array4d &TD, const double TR, const int N, const Eigen::Array3d flip);
        MP3RAGE(std::istream &istr, const bool prompt);
        size_t size() const override { return 3; }
        Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override { QI_EXCEPTION("Not implemented"); }
        Eigen::ArrayXcd signal(const double M0, const double T1, const double B1, const double eta) const;
        void write(std::ostream &os) const override;
        std::string name() const override { return "MP3RAGE"; }
};

} // End namespace QI

#endif // SEQUENCES_MPRAGE_H
