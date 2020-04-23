// #define QI_DEBUG_BUILD 1
#include "Macro.h"
#include "rufis_ss.hpp"
#include "ss_mt.h"

using AugMat = Eigen::Matrix<double, 5, 5>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 5>;

auto SS_MT_Model::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T = double;

    T const &M0_f = v[0];
    T const &M0_b = v[1];
    T const &R1_f = 1. / v[2];
    T const &R1_b = R1_f;
    T const &R2_f = 1. / v[3];
    T const &T2_b = v[4];
    T const &k    = v[5];
    T const &k_bf = k * M0_f / (M0_f + M0_b);
    T const &k_fb = k * M0_b / (M0_f + M0_b);
    T const &f0   = v[6];
    T const &B1p  = v[7];

    AugMat R;
    R << -R2_f, 0, 0, 0, 0,          //
        0, -R2_f, 0, 0, 0,           //
        0, 0, -R1_f, 0, M0_f * R1_f, //
        0, 0, 0, -R1_b, M0_b * R1_b, //
        0, 0, 0, 0, 0;

    AugMat K;
    K << 0, 0, 0, 0, 0,       //
        0, 0, 0, 0, 0,        //
        0, 0, -k_fb, k_bf, 0, //
        0, 0, k_fb, -k_bf, 0, //
        0, 0, 0, 0, 0;

    AugMat const S = Eigen::DiagonalMatrix<double, 5, 5>({0, 0, 1., 1., 1.}).toDenseMatrix();

    AugMat const RpK = R + K;

    // Setup constant matrices
    AugMat const Rrd  = (RpK * (sequence.TR - sequence.Trf)).exp();
    AugMat const ramp = (RpK * sequence.Tramp).exp();

    auto RF = [&RpK, &T2_b, &f0, &B1p, this](double const alpha,
                                             double const tau,
                                             double const df,
                                             double const p1,
                                             double const p2) {
        T const B1 = B1p * alpha / (p1 * tau);
        T const dw = 2. * M_PI * (f0 + df);

        T const G = this->lineshape(f0 + df, T2_b);
        T const W = M_PI * B1p * B1p * G * (p2 / (p1 * p1)) * (alpha * alpha) / (tau * tau);

        AugMat rf;
        rf << 0, dw, 0, 0, 0, //
            -dw, 0, B1, 0, 0, //
            0, -B1, 0, 0, 0,  //
            0, 0, 0, -W, 0,   //
            0, 0, 0, 0, 0;
        QI_DBMAT(rf);
        AugMat const Arf = ((rf + RpK) * tau).exp();
        QI_DBMAT(Arf);
        return Arf;
    };

    Eigen::ArrayXd sig(sequence.size());
    for (long is = 0; is < sequence.size(); is++) {
        AugMat const rfp     = RF(sequence.prep_FA[is],
                              sequence.prep_Trf,
                              sequence.prep_df[is],
                              sequence.prep_p1,
                              sequence.prep_p2);
        AugMat const rf1     = RF(sequence.FA[is], sequence.Trf, 0., 1., 1.);
        AugMat       TR_mat  = S * Rrd * rf1;
        AugMat       seg_mat = TR_mat.pow(sequence.spokes_per_seg);

        // Calculate the steady-state just before the segment readout
        AugMat X    = AugMat::Identity();
        X           = ramp * S * rfp * ramp * seg_mat * X;
        AugVec m_ss = SolveSteadyState(X);

        auto m_gm = GeometricAvg(TR_mat, seg_mat, m_ss, sequence.spokes_per_seg);
        sig[is]   = m_gm[2] * sin(B1p * sequence.FA[is]);
    }
    QI_DBVEC(sig)
    return sig;
}

void SS_MT_Model::derived(const VaryingArray &varying,
                          const FixedArray & /* Unused */,
                          DerivedArray &derived) const {
    const auto &M0_f = varying[0];
    const auto &M0_b = varying[1];
    derived[0]       = 100.0 * M0_b / (M0_f + M0_b);
}