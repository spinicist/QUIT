// #define QI_DEBUG_BUILD 1
#include "transient_qmt_model.h"
#include "Macro.h"
#include "parmesan.hpp"

using AugMat = Eigen::Matrix<double, 6, 6>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 6>;

auto QMTModel::signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(double) {
    using T = double;

    T const M0   = v[0];
    T const R1_f = 1. / v[1];
    T const f_s  = v[2];
    T const f_f  = 1.0 - f_s;
    T const B1   = v[3];
    T const R2_f = 1. / f[0];
    T const R1_s = 1. / f[1];
    T const R2_s = 1. / f[2];
    T const Rx   = f[3];

    T const M0_s = f_s * M0;
    T const M0_f = M0 - M0_s;

    T const R_fs = Rx * f_s;
    T const R_sf = Rx * f_f;

    QI_DBVEC(v)
    QI_DB(M0)
    QI_DB(f_f)
    QI_DB(f_s)
    QI_DB(M0_f)
    QI_DB(M0_s)
    QI_DB(R1_f)
    QI_DB(R2_f)
    QI_DB(R1_s)
    QI_DB(R2_s)
    QI_DB(Rx)
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
    QI_DBMAT(R);
    QI_DBMAT(K);
    QI_DBMAT(RpK);
    using DT       = Eigen::DiagonalMatrix<double, 6, 6>;
    AugMat const S = DT(DT::DiagonalVectorType{0, 0, 1., 0., 1., 1.}).toDenseMatrix();
    QI_DBMAT(S);
    // Setup readout segment matrices
    AugMat const ramp = S * (RpK * sequence.Tramp).exp();
    QI_DBMAT(ramp);
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    for (int ip = 0; ip < sequence.FA.rows(); ip++) {
        double const B1x  = B1 * sequence.FA[ip] / sequence.Trf;
        double const r2sl = 0;
        AugMat       rf;
        rf << 0, 0, 0, 0, 0, 0,                          //
            0, 0, B1x, 0, 0, 0,                          //
            0, -B1x, 0, 0, 0, 0,                         //
            0, 0, 0, -r2sl, B1x, 0, 0, 0, 0, -B1x, 0, 0, //
            0, 0, 0, 0, 0, 0;
        A_mats[ip]   = ((RpK + rf) * sequence.Trf).exp();
        R_mats[ip]   = (RpK * (sequence.TR - sequence.Trf)).exp();
        seg_mats[ip] = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS);
    }

    // Setup pulse matrices
    std::vector<AugMat> prep_mats(sequence.preps());
    for (int ip = 0; ip < sequence.preps(); ip++) {
        double const B1x  = B1 * sequence.FAprep[ip] / sequence.Tprep;
        double const r2sl = 8e3;
        AugMat       rf;
        rf << 0, 0, 0, 0, 0, 0,     //
            0, 0, B1x, 0, 0, 0,     //
            0, -B1x, 0, 0, 0, 0,    //
            0, 0, 0, -r2sl, B1x, 0, //
            0, 0, 0, -B1x, 0, 0,    //
            0, 0, 0, 0, 0, 0;
        QI_DBMAT(RpK);
        QI_DBMAT(rf);
        QI_DB(sequence.Tprep);
        prep_mats[ip] = ((RpK + rf) * sequence.Tprep).exp();
        QI_DBMAT(prep_mats[ip]);
    }

    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int ip = 0; ip < sequence.preps(); ip++) {
        X = ramp * seg_mats[ip] * ramp * prep_mats[ip] * X;
    }
    AugVec const m_ss = SolveSteadyState(X);

    // Now loop through the segments and record the signal for each
    Eigen::ArrayXd sig(sequence.SPS * sequence.FAprep.size());
    QI_DBVEC(m_ss);
    AugVec m  = m_ss;
    int    ii = 0;
    for (int ip = 0; ip < sequence.preps(); ip++) {
        QI_DBVEC(m);
        QI_DBVEC((prep_mats[ip] * m))
        m = ramp * prep_mats[ip] * m;
        QI_DBVEC(m);
        for (int is = 0; is < sequence.SPS; is++) {
            m = A_mats[ip] * m;
            QI_DBVEC(m);
            sig[ii++] = m[1];
            m         = S * R_mats[ip] * m;
        }
        m = ramp * m;
    }
    QI_DBVEC(sig);
    if (sequence.basis.size()) {
        return sequence.basis * sig.matrix();
    } else {
        return sig;
    }
}
