// #define QI_DEBUG_BUILD 1
#include "Macro.h"
#include "mt_model.h"
#include "rufis_ss.hpp"

using AugMat = Eigen::Matrix<double, 5, 5>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 5>;

auto MTModel::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
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
    T const &B1   = v[6];

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
    AugMat const Rrd          = (RpK * (sequence.TR - sequence.Trf)).exp();
    AugMat const ramp         = (RpK * sequence.Tramp).exp();
    AugMat const prep_spoiler = S * (RpK * sequence.Tspoil).exp();

    Eigen::ArrayXd sig(sequence.size());
    for (long is = 0; is < sequence.RUFIS_FA.rows(); is++) {
        double const rf1_B1x = B1 * sequence.RUFIS_FA[is] / sequence.Trf;
        double const rf1_W   = M_PI * lineshape({0}, T2_b) * rf1_B1x * rf1_B1x;
        AugMat       rf1;
        rf1 << 0, 0, 0, 0, 0,     //
            0, 0, rf1_B1x, 0, 0,  //
            0, -rf1_B1x, 0, 0, 0, //
            0, 0, 0, -rf1_W, 0,   //
            0, 0, 0, 0, 0;

        AugMat const Ard     = ((RpK + rf1) * sequence.Trf).exp();
        AugMat       TR_mat  = S * Rrd * Ard;
        AugMat       seg_mat = TR_mat.pow(sequence.SPS);

        auto const & p        = sequence.MT_pulse;
        AugMat       MT       = AugMat::Identity();
        double const dt_act   = sequence.MT_pulsewidth / p.B1x.size();
        double const MT_scale = (sequence.MT_FA[is] / p.FA) * (p.width / sequence.MT_pulsewidth);
        double const dw       = 2 * M_PI * sequence.MT_offsets[is];
        for (long im = 0; im < p.B1x.size(); im++) {
            double const gamma  = 267.52219; // radians per second per uT
            double const MT_B1x = gamma * B1 * p.B1x[im] * MT_scale;
            double const MT_W   = M_PI * lineshape(sequence.MT_offsets[is], T2_b) * MT_B1x * MT_B1x;
            AugMat       MT_rf;
            MT_rf << 0, dw, 0, 0, 0,  //
                -dw, 0, MT_B1x, 0, 0, //
                0, -MT_B1x, 0, 0, 0,  //
                0, 0, 0, -MT_W, 0,    //
                0, 0, 0, 0, 0;
            AugMat const MT_step = ((RpK + MT_rf) * dt_act).exp();
            MT                   = MT_step * MT;
        }
        QI_DB(dt_act)
        QI_DB(MT_scale)
        QI_DB(sequence.MT_FA[is])
        QI_DB(p.FA)
        QI_DB(p.width)
        QI_DB(sequence.MT_pulsewidth)
        QI_DB(sequence.MT_FA[is])
        QI_DB(sequence.MT_offsets[is])
        QI_DBMAT(MT)
        // Calculate the steady-state just before the segment readout
        AugMat X    = AugMat::Identity();
        X           = ramp * prep_spoiler * MT * ramp * seg_mat * X;
        AugVec m_ss = SolveSteadyState(X);

        auto       m_gm   = GeometricAvg(TR_mat, seg_mat, m_ss, sequence.SPS);
        auto const signal = m_gm[2] * sin(B1 * sequence.RUFIS_FA[is]);
        sig[is]           = signal;
    }
    QI_DBVEC(sig)
    return sig;
}

void MTModel::derived(const VaryingArray &varying,
                      const FixedArray & /* Unused */,
                      DerivedArray &derived) const {
    const auto &M0_f = varying[0];
    const auto &M0_b = varying[1];
    derived[0]       = 100.0 * M0_b / (M0_f + M0_b);
}