// #define QI_DEBUG_BUILD 1
#include "mupa_model_b1.h"
#include "Macro.h"
#include "rufis_ss.hpp"

using AugMat = Eigen::Matrix<double, 4, 4>;
using AugVec = Eigen::Vector<double, 4>;

auto MUPAB1Model::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T = double;

    T const &M0 = v[0];
    T const &R1 = 1. / v[1];
    T const &R2 = 1. / v[2];
    T const &B1 = v[3];

    QI_DBVEC(v);
    QI_DB(M0);
    QI_DB(R1);
    QI_DB(R2);
    QI_DB(B1);
    QI_DBVEC(sequence.Trf)
    QI_DBVEC(sequence.FA)
    QI_DB(sequence.SPS)

    AugMat R;
    R << -R2, 0, 0, 0,      //
        0, -R2, 0, 0,       //
        0, 0, -R1, M0 * R1, //
        0, 0, 0, 0;

    AugMat const Rrd  = (R * sequence.TR).exp();
    AugMat const S    = Eigen::DiagonalMatrix<double, 4, 4>({0, 0, 1., 1.}).toDenseMatrix();
    AugMat const ramp = (R * sequence.Tramp).exp();

    // Setup readout segment matrices
    std::vector<AugMat> TR_mats(sequence.size());
    std::vector<AugMat> seg_mats(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        AugMat rf;
        float  B1x = B1 * sequence.FA[is] / sequence.Trf[is];
        rf << 0, 0, 0, 0,  //
            0, 0, B1x, 0,  //
            0, -B1x, 0, 0, //
            0, 0, 0, 0;
        AugMat const Ard = ((R + rf) * sequence.Trf[is]).exp();
        TR_mats[is]      = S * Rrd * Ard;
        seg_mats[is]     = TR_mats[is].pow(sequence.SPS[is]);
    }

    // Setup pulse matrices
    std::vector<AugMat> prep_mats(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        auto const & name = sequence.prep[is];
        auto const & p    = sequence.prep_pulses[name];
        double const E2   = exp(-R2 * p.T_trans);
        double const E1   = exp(-R1 * p.T_long);
        AugMat       C;
        C << 0, 0, 0, 0,                                 //
            0, 0, 0, 0,                                  //
            0, 0, E1 * E2 * cos(p.FAeff), M0 * (1 - E1), //
            0, 0, 0, 1;
        prep_mats[is] = C;
    }

    // First calculate the system matrix and get SS
    QI_DBMAT(R);
    AugMat X = AugMat::Identity();
    for (int is = 0; is < sequence.size(); is++) {
        X = ramp * seg_mats[is] * ramp * prep_mats[is] * X;
    }
    AugVec m_aug = SolveSteadyState(X);
    QI_DBMAT(X);
    QI_DBVEC(m_aug);
    // Now loop through the segments and record the signal for each
    Eigen::ArrayXd sig(sequence.size());
    for (int is = 0; is < sequence.size(); is++) {
        m_aug     = ramp * prep_mats[is] * m_aug;
        auto m_gm = GeometricAvg(TR_mats[is], seg_mats[is], m_aug, sequence.SPS[is]);
        sig[is]   = m_gm[2] * sin(B1 * sequence.FA[is]);
        m_aug     = ramp * seg_mats[is] * m_aug;
    }
    QI_DBVEC(sig);
    return sig;
}