// #define QI_DEBUG_BUILD 1
#include "transient_qmt_model.h"
#include "Macro.h"
#include "parmesan.hpp"

using AugMat = Eigen::Matrix<double, 6, 6>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 6>;

auto QMTModel::signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(double) {
    using T = double;

    T const &M0   = v[0];
    T const &R1_f = 1. / v[1];
    T const &f_s  = v[2];
    T const &B1   = v[3];
    T const &R2_f = 1. / f[0];
    T const &R1_s = 1. / f[1];
    T const &T2_s = f[2];
    T const &k    = f[3];

    T const M0_s = f_s * M0;
    T const M0_f = M0 - M0_s;

    T const &R_sf = k / f_s;
    T const &R_fs = k / (1 - f_s);

    QI_DBVEC(v)
    QI_DB(M0_f)
    QI_DB(M0_s)
    QI_DB(R1_f)
    QI_DB(R2_f)
    QI_DB(R1_s)
    QI_DB(R2_s)
    QI_DB(R_sf)
    QI_DB(R_fs)
    QI_DB(B1)

    // State vector is [x_f y_f z_f x_s z_s 1]

    AugMat R;
    R << -R2_f, 0, 0, 0, 0, 0,          //
        0, -R2_f, 0, 0, 0, 0,           //
        0, 0, -R1_f, 0, 0, M0_f * R1_f, //
        0, 0, 0, 0, 0, 0,               //
        0, 0, 0, 0, -R1_s, M0_s * R1_s, //
        0, 0, 0, 0, 0, 0;
    AugMat K;
    K << 0, 0, 0, 0, 0, 0,       //
        0, 0, 0, 0, 0, 0,        //
        0, 0, -R_fs, 0, R_sf, 0, //
        0, 0, 0, 0, 0, 0,        //
        0, 0, R_fs, 0, -R_sf, 0, //
        0, 0, 0, 0, 0, 0;
    AugMat const RpK = R + K;
    using DT         = Eigen::DiagonalMatrix<double, 6, 6>;
    AugMat const S   = DT(DT::DiagonalVectorType{0, 0, 1., 0., 1., 1.}).toDenseMatrix();

    // Setup readout segment matrices
    AugMat const        ramp = (RpK * sequence.Tramp).exp();
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    for (int is = 0; is < sequence.FA.rows(); is++) {
        double const B1x  = B1 * sequence.FA[is] / sequence.Trf;
        double const r2sl = 0;
        AugMat       rf;
        rf << 0, 0, 0, 0, 0, 0,                          //
            0, 0, B1x, 0, 0, 0,                          //
            0, -B1x, 0, 0, 0, 0,                         //
            0, 0, 0, -r2sl, B1x, 0, 0, 0, 0, -B1x, 0, 0, //
            0, 0, 0, 0, 0, 0;
        A_mats[is]   = ((RpK + rf) * sequence.Trf).exp();
        R_mats[is]   = (RpK * (sequence.TR - sequence.Trf)).exp();
        seg_mats[is] = (R_mats[is] * A_mats[is]).pow(sequence.SPS);
    }

    // Setup pulse matrices
    std::vector<AugMat> prep_mats(sequence.preps());
    for (int ip = 0; ip < sequence.preps(); ip++) {
        double const B1x  = B1 * sequence.FAprep[ip] / sequence.Tprep;
        double const r2sl = 8e3;
        AugMat       rf;
        rf << 0, 0, 0, 0, 0, 0,                          //
            0, 0, B1x, 0, 0, 0,                          //
            0, -B1x, 0, 0, 0, 0,                         //
            0, 0, 0, -r2sl, B1x, 0, 0, 0, 0, -B1x, 0, 0, //
            0, 0, 0, 0, 0, 0;
        prep_mats[ip] = ((RpK + rf) * sequence.Tprep).exp();
    }

    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int ip = 0; ip < sequence.preps(); ip++) {
        X = ramp * seg_mats[ip] * ramp * S * prep_mats[ip] * X;
    }
    AugVec m_ss = SolveSteadyState(X);

    // Now loop through the segments and record the signal for each
    Eigen::ArrayXd sig(sequence.size());
    QI_DBVEC(m_ss);
    AugVec m  = m_ss;
    int    ii = 0;
    for (int ip = 0; ip < sequence.preps(); ip++) {
        m = ramp * S * prep_mats[ip] * m;
        for (int is = 0; is < sequence.SPS; is++) {
            m         = A_mats[ip] * m;
            sig[ii++] = m[1];
            m         = S * R_mats[ip] * m;
        }
        m = ramp * m;
    }
    QI_DBVEC(v);
    if (sequence.basis.size()) {
        return sequence.basis * sig.matrix();
    } else {
        return sig;
    }
}
