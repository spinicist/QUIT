// #define QI_DEBUG_BUILD 1
#include "Macro.h"
#include "parmesan.hpp"
#include "ss_T2.h"

auto SS_T1T2_Model::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T      = double;
    using AugMat = Eigen::Matrix<T, 4, 4>; // Short for Augmented Matrix
    using AugVec = Eigen::Vector<T, 4>;

    QI_DBVEC(v);
    T const &M0     = v[0];
    T const &R1     = 1. / v[1];
    T const &R2     = 1. / v[2];
    T const &f0     = v[3];
    T const &B1plus = v[4];

    // Relaxation
    AugMat R;
    R << -R2, 0, 0, 0, //
        0, -R2, 0, 0,  //
        0, 0, -R1, R1, //
        0, 0, 0, 0;

    // Spoiling
    AugMat S;
    S << 0, 0, 0, 0, //
        0, 0, 0, 0,  //
        0, 0, 1, 0,  //
        0, 0, 0, 1;

    QI_DBMAT(R);
    // Useful for later
    auto RF =
        [&R, &f0, &B1plus](double const alpha, double const tau, double const df, double const p1) {
            AugMat  rf;
            T const B1nom = alpha / (p1 * tau);
            T const B1    = B1plus * B1nom;
            T const dw    = 2. * M_PI * (f0 + df);
            rf << 0, dw, 0, 0, //
                -dw, 0, B1, 0, //
                0, -B1, 0, 0,  //
                0, 0, 0, 0;
            QI_DBMAT(rf);
            AugMat const Arf = ((rf + R) * tau).exp();
            QI_DBMAT(Arf);
            return Arf;
        };

    // Setup constant matrices
    AugMat const Rrd  = (R * (sequence.TR - sequence.Trf)).exp();
    AugMat const ramp = (R * sequence.Tramp).exp();
    // AugMat const prep_spoiler = (R * sequence.Tspoil).exp();

    Eigen::ArrayXd sig(sequence.size());
    for (long is = 0, ie = sequence.size(); is < ie; is++) {
        AugMat const rfp =
            RF(sequence.prep_FA[is], sequence.prep_Trf, sequence.prep_df[is], sequence.prep_p1);
        AugMat const rf1     = RF(sequence.FA[is], sequence.Trf, 0., 1.);
        AugMat const TR_mat  = S * Rrd * rf1;
        AugMat const seg_mat = TR_mat.pow(sequence.spokes_per_seg);

        // Calculate the steady-state just before the segment readout
        AugMat const X    = ramp * S * rfp * ramp * seg_mat;
        AugVec const m_ss = SolveSteadyState(X);

        auto       m_gm   = GeometricAvg(TR_mat, seg_mat, m_ss, sequence.spokes_per_seg);
        auto const signal = m_gm[2] * sin(B1plus * sequence.FA[is]);
        sig[is]           = M0 * signal;
    }
    QI_DBVEC(v)
    QI_DBVEC(sig)
    return sig;
}