// #define QI_DEBUG_BUILD 1
#include "ss_model.h"
#include "Macro.h"
#include "rufis_ss.hpp"

auto SS_T1_Model::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T      = double;
    using AugMat = Eigen::Matrix<T, 2, 2>; // Short for Augmented Matrix
    using AugVec = Eigen::Vector<T, 2>;

    T const &M0 = v[0];
    T const &R1 = 1. / v[1];
    T const &B1 = v[2];

    AugMat R;
    R << -R1, M0 * R1, //
        0, 0;

    // Setup constant matrices
    AugMat const Rrd  = (R * sequence.TR).exp();
    AugMat const ramp = (R * sequence.Tramp).exp();
    // AugMat const prep_spoiler = (R * sequence.Tspoil).exp();

    Eigen::ArrayXd sig(sequence.size());
    for (long is = 0; is < sequence.size(); is++) {
        AugMat rf1;
        rf1 << cos(B1 * sequence.FA[is]), 0, //
            0, 1;

        AugMat TR_mat  = Rrd * rf1;
        AugMat seg_mat = TR_mat.pow(sequence.spokes_per_seg);

        AugMat rfp;
        rfp << cos(B1 * sequence.prep_FA[is]), 0, //
            0, 1;

        // Calculate the steady-state just before the segment readout
        AugMat X    = AugMat::Identity();
        X           = ramp * rfp * ramp * seg_mat * X;
        AugVec m_ss = SolveSteadyState(X);

        auto       m_gm   = GeometricAvg(TR_mat, seg_mat, m_ss, sequence.spokes_per_seg);
        auto const signal = m_gm[0] * sin(B1 * sequence.FA[is]);
        sig[is]           = signal;
        QI_DBMAT(rf1)
        QI_DBMAT(rfp)
        QI_DBVEC(m_ss)
        QI_DBVEC(m_gm)
    }
    QI_DBVEC(v)
    QI_DBVEC(sig)
    return sig;
}
