/*
 *  MPRAGESequence.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MPRAGESequence.h"
#include "IO.h"

using namespace std;
using namespace Eigen;

namespace QI {

MPRAGE::MPRAGE(const ArrayXd &TI, const ArrayXd &TD, const double TR, const int Nseg, const int Nk0, const double flip, const double eta) :
    SequenceBase(), m_TI(TI), m_TD(TD), m_ETL(Nseg), m_k0(Nk0), m_eta(eta) {
    m_TR = TR;
    m_flip.resize(1); m_flip[0] = flip;
    if (m_TI.size() != m_TD.size()) {
        QI_EXCEPTION("Inversion and delay time arrays must match in size.");
    }
}

MPRAGE::MPRAGE(std::istream &istr, const bool prompt) : SequenceBase() {
    if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
    double inFlip, TRseg;
    QI::Read(istr, inFlip);
    m_flip = ArrayXd::Ones(1) * inFlip * M_PI / 180.;
    if (prompt) cout << "Enter read-out TR (seconds): " << flush;
    QI::Read(istr, m_TR);
    if (prompt) cout << "Enter segment size: " << flush;
    QI::Read(istr, m_ETL);
    if (prompt) cout << "Enter k0: " << flush;
    QI::Read(istr, m_k0);
    if (prompt) cout << "Enter inversion times (seconds): " << flush;
    QI::ReadArray(istr, m_TI);
    if (prompt) cout << "Enter relaxation delay times (seconds): " << flush;
    QI::ReadArray(istr, m_TD);
    if (m_TI.rows() != m_TD.rows()) {
        QI_EXCEPTION("Must have the same number of relaxation delays as inversion times.");
    }
    if (prompt) cout << "Enter inversion efficiency (<= 1.0): " << flush;
    QI::Read(istr, m_eta);
    if (m_eta > 1.0) {
        QI_EXCEPTION("Inversion-efficiency cannot exceed 1.0");
    }
}

IRSPGR::IRSPGR(std::istream &istr, const bool prompt) : MPRAGE() {
    if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
    double inFlip;
    QI::Read(istr, inFlip);
    m_flip = ArrayXd::Ones(1) * inFlip * M_PI / 180.;
    if (prompt) cout << "Enter read-out TR (seconds): " << flush;
    QI::Read(istr, m_TR);

    int NPE2;
    if (prompt) cout << "Enter number of spatial locations (remember +4): ";
    QI::Read(istr, NPE2);
    if (NPE2 >= 64) {
     m_ETL = NPE2 / 2;
    } else {
     m_ETL = NPE2;
    }
    m_k0 = 0;
    
    if (prompt) cout << "Enter TIs (seconds): " << flush;
    QI::ReadArray(istr, m_TI);

    m_TD = ArrayXd::Zero(m_TI.size()); // For GE IR-SPGR the delay time is zero
    m_eta = 1.0; // Assume this for now
}

ArrayXcd MPRAGE::signal(shared_ptr<Model> m, const VectorXd &par) const {
    return m->MPRAGE(par, m_flip[0], m_TR, m_ETL, m_k0, m_eta, m_TI, m_TD);
}

void MPRAGE::write(ostream &os) const {
    os << name() << endl;
    os << "TR: " << m_TR << "\tSegment Length: " << m_ETL << "\tk-Zero: " << m_k0 << "\tAlpha: " << m_flip[0] * 180 / M_PI << endl;
    os << "TI: " << m_TI.transpose() << "\tTD: " << m_TD.transpose() << "\tEta: " << m_eta << "\tSegment Time: " << (m_TI + m_TR*(m_ETL - m_k0) + m_TD) << endl;
}

ArrayXd MPRAGE::weights(const double f0) const {
    return ArrayXd::Ones(size()) * 1.0;
}

MP2RAGE::MP2RAGE(const Array3d &TD, const double TR, const int N, const Array2d flip) :
    SequenceBase(), m_TD(TD), m_N(N) {
    m_TR = TR;
    m_flip = flip;
}

MP2RAGE::MP2RAGE(std::istream &istr, const bool prompt) : SequenceBase() {
    if (prompt) cout << "Enter read-out flip-angles (degrees): " << flush;
    QI::ReadArray(istr, m_flip);
    if (m_flip.size() != 2) {
        QI_EXCEPTION("Must have 2 flip-angles for MP3RAGE");
    }
    m_flip *= M_PI / 180.;
    if (prompt) cout << "Enter read-out TR (seconds): " << flush;
    QI::Read(istr, m_TR);
    if (prompt) cout << "Enter segment size: " << flush;
    QI::Read(istr, m_N);
    if (prompt) cout << "Enter inversion times (seconds): " << flush;
    ArrayXd m_TI;
    QI::ReadArray(istr, m_TI);
    if (m_TI.size() != 2) {
        QI_EXCEPTION("Must specify 2 TI times for MP3RAGE");
    }
    // Assume centric for now
    m_TD[0] = m_TI[0];
    m_TD[1] = m_TI[1] - m_TD[0] - m_N*m_TR;
    if (prompt) cout << "Enter overall TR (seconds): " << flush;
    float TRseg;
    QI::Read(istr, TRseg);
    m_TD[2] = TRseg - m_TD[1] - m_N*m_TR;
}

ArrayXcd MP2RAGE::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP2RAGE(m_flip, m_TR, m_N, m_TD, M0, T1, B1, eta);
}

void MP2RAGE::write(ostream &os) const {
    os << name() << endl;
    os << "TR: " << m_TR << "\tN: " << m_N << "\tAlpha: " << m_flip[0] * 180 / M_PI << endl;
    os << "TD: " << m_TD.transpose() << "\tTS: " << (m_TD.sum() + 2*m_N*m_TR) << endl;
}

MP3RAGE::MP3RAGE(const Array4d &TD, const double TR, const int N, const Array3d flip) :
    SequenceBase(), m_TD(TD), m_N(N) {
    m_TR = TR;
    m_flip = flip;
}

MP3RAGE::MP3RAGE(std::istream &istr, const bool prompt) : SequenceBase() {
    if (prompt) cout << "Enter read-out flip-angles (degrees): " << flush;
    QI::ReadArray(istr, m_flip);
    if (m_flip.size() != 3) {
        QI_EXCEPTION("Must have 3 flip-angles for MP3RAGE");
    }
    m_flip *= M_PI / 180.;
    if (prompt) cout << "Enter read-out TR (seconds): " << flush;
    QI::Read(istr, m_TR);
    if (prompt) cout << "Enter segment size: " << flush;
    QI::Read(istr, m_N);
    if (prompt) cout << "Enter inversion times (seconds): " << flush;
    ArrayXd m_TI;
    QI::ReadArray(istr, m_TI);
    if (m_TI.size() != 3) {
        QI_EXCEPTION("Must specify 3 TI times for MP3RAGE");
    }
    // Assume centric for now
    m_TD[0] = m_TI[0];
    m_TD[1] = m_TI[1] - m_TD[0] - m_N*m_TR;
    m_TD[2] = m_TI[2] - m_TD[1] - m_N*m_TR;
    if (prompt) cout << "Enter overall TR (seconds): " << flush;
    float TRseg;
    QI::Read(istr, TRseg);
    m_TD[3] = TRseg - m_TD[2] - m_N*m_TR;
}

ArrayXcd MP3RAGE::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP3RAGE(m_flip, m_TR, m_N, m_TD, M0, T1, B1, eta);
}

void MP3RAGE::write(ostream &os) const {
    os << name() << endl;
    os << "TR: " << m_TR << "\tN: " << m_N << "\tAlpha: " << m_flip[0] * 180 / M_PI << endl;
    os << "TD: " << m_TD.transpose() << "\tTS: " << (m_TD.sum() + 3*m_N*m_TR) << endl;
}

} // End namespace QI