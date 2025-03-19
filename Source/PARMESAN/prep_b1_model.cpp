// #define QI_DEBUG_BUILD 1
#include "prep_b1_model.h"
#include "Macro.h"
#include "parmesan.hpp"

using AugMat = Eigen::Matrix<double, 3, 3>;
using AugVec = Eigen::Vector<double, 3>;

auto PrepB1Model::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
    using T = double;

    T const      M0 = v[0];
    T const      R1 = 1. / v[1];
    double const R2 = 10.;
    T const      B1 = v[2];

    QI_DBVEC(v);
    QI_DB(M0);
    QI_DB(R1);
    QI_DB(R2);
    QI_DB(B1);
    QI_DBVEC(sequence.Trf)
    QI_DBVEC(sequence.FA)
    QI_DB(sequence.spokes_per_seg)
    QI_DBVEC(sequence.groups_per_seg)

    AugMat R;
    R << -R2, 0, 0, //
        0, -R1, R1, //
        0, 0, 0;

    AugMat const Rrd  = (R * sequence.TR).exp();
    using DT          = Eigen::DiagonalMatrix<double, 3, 3>;
    AugMat const S    = DT(DT::DiagonalVectorType{0., 1., 1.}).toDenseMatrix();
    AugMat const ramp = S * (R * sequence.Tramp).exp();

    // Setup pulse matrices
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    std::vector<AugMat> prep_mats(sequence.preps());
    for (int ip = 0; ip < sequence.preps(); ip++) {
        double const B1x = B1 * sequence.FA[ip] / sequence.Trf;
        AugMat       rf;
        rf << 0, B1x, 0, //
            -B1x, 0, 0,  //
            0, 0, 0;
        A_mats[ip]   = ((R + rf) * sequence.Trf).exp();
        R_mats[ip]   = (R * (sequence.TR - sequence.Trf)).exp();
        seg_mats[ip] = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS);

        double const B1p = B1 * sequence.FAprep[ip] / sequence.Tprep;
        AugMat       rfp;
        rfp << 0, B1p, 0, //
            -B1p, 0, 0,   //
            0, 0, 0;
        QI_DBMAT(RpK);
        QI_DBMAT(rf);
        QI_DB(sequence.Tprep);
        prep_mats[ip] = ((R + rfp) * sequence.Tprep).exp();
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