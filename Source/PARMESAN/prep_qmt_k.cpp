#include "prep_qmt_k.h"

namespace {
using AugMat = Eigen::Matrix<double, 6, 6>; // Short for Augmented Matrix
using AugVec = Eigen::Vector<double, 6>;

auto Relax(double const M0_f,
           double const R1_f,
           double const R2_f,
           double const M0_s,
           double const R1_s,
           double const R2_s) -> AugMat {
    AugMat R;
    R.setZero();
    R(0, 0) = -R2_f;
    R(1, 1) = -R2_f;
    R(2, 2) = -R1_f;
    R(2, 5) = M0_f * R1_f;
    R(3, 3) = -R2_s;
    R(4, 4) = -R1_s;
    R(4, 5) = M0_s * R1_s;
    return R;
}

auto Exchange(double const R_fs, double const R_sf) -> AugMat {
    AugMat K;
    K.setZero();
    K(2, 2) = -R_fs;
    K(2, 4) = R_sf;
    K(4, 2) = R_fs;
    K(4, 4) = -R_sf;
    return K;
}

auto RF(double const       flip,
        double const       pw,
        double const       B1,
        double const       R2_s,
        QI::RegularGrid const &A_sl) -> AugMat {
    double const B1x = B1 * flip / pw;
    double const τ   = pw * R2_s;
    double const ω   = B1 * flip / τ;
    AugMat       rf;
    rf.setZero();
    rf(1, 2) = B1x;
    rf(2, 1) = -B1x;
    rf(3, 3) = -R2_s * A_sl(ω, τ);
    rf(3, 4) = B1x;
    rf(4, 3) = -B1x;
    return rf;
}

auto Spoil() -> AugMat {
    AugMat S;
    S.setZero();
    S(2, 2) = 1;
    S(4, 4) = 1;
    S(5, 5) = 1;
    return S;
}
} // namespace

namespace QI {

int PrepQMTk::input_size(const int /* Unused */) const {
    return sequence.size();
}

auto PrepQMTk::signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType) {
    double const M0   = v[0];
    double const R1_f = 1. / v[1];
    double const f_s  = v[2];
    double const f_f  = 1.0 - f_s;
    double const B1   = v[3];
    double const R2_f = 1. / f[0];
    double const R1_s = 1. / f[1];
    double const R2_s = 1. / f[2];
    double const R_fs = f[3] / f_f;
    double const R_sf = f[3] / f_s;

    // State vector is [x_f y_f z_f x_s z_s 1]
    AugMat       R        = Relax(f_f, R1_f, R2_f, f_s, R1_s, R2_s);
    AugMat       K        = Exchange(R_fs, R_sf);
    AugMat const RpK      = R + K;
    AugMat const S        = Spoil();
    AugMat const ramp     = S * (RpK * sequence.Tramp).exp();
    AugMat const spoilTRs = S * (RpK * sequence.TR * sequence.spoilers).exp();
    AugMat const Dprep    = (RpK * sequence.Dprep).exp();
    AugMat const Dseg     = (RpK * sequence.Dseg).exp();
    // QI_DBMAT(R)
    // QI_DBMAT(K)
    // QI_DBMAT(S)
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    std::vector<AugMat> prep_mats(sequence.preps());
    for (int ip = 0; ip < sequence.FA.rows(); ip++) {
        auto const rfa = RF(sequence.FA[ip], sequence.Trf, B1, R2_s, A_sl);
        A_mats[ip]     = ((RpK + rfa) * sequence.Trf).exp();
        R_mats[ip]     = (RpK * (sequence.TR - sequence.Trf)).exp();
        seg_mats[ip]   = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS) * spoilTRs;
        AugMat rfp     = RF(sequence.FAprep[ip], sequence.Tprep[ip], B1, R2_s, A_sl);
        prep_mats[ip]  = ((RpK + rfp) * sequence.Tprep[ip]).exp();
    }
    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int ip = 0; ip < sequence.preps(); ip++) {
        X = Dseg * ramp * seg_mats[ip] * ramp * Dprep * prep_mats[ip] * X;
    }
    AugVec const m_ss = SolveSteadyState(X);
    QI_DBVEC(m_ss);
    // Now loop through the segments and record the signal for each
    QI_ARRAY(DataType) sig(sequence.SPS * sequence.FAprep.size());
    AugVec m  = m_ss;
    int    ii = 0;
    for (int ip = 0; ip < sequence.preps(); ip++) {
        m = spoilTRs * ramp * Dprep * prep_mats[ip] * m;
        for (int is = 0; is < sequence.SPS; is++) {
            m         = A_mats[ip] * m;
            sig[ii++] = M0 * m[1];
            m         = S * R_mats[ip] * m;
        }
        m = Dseg * ramp * m;
    }
    if (basis.size()) {
        return basis * sig.matrix();
    } else {
        return sig;
    }
}

int PrepQMTkFull::input_size(const int /* Unused */) const {
    return sequence.size();
}

auto PrepQMTkFull::signal(VaryingArray const &v, FixedArray const &f) const
    -> QI_ARRAY(DataType) {
    // "M0", "T1_f", "T2_f", "T1_s", "T2_s", "f_s", "Rx", "B1"
    double const M0   = v[0];
    double const R1_f = 1. / v[1];
    double const R2_f = 1. / f[0];
    double const R1_s = 1. / v[2];
    double const R2_s = 1. / v[3];
    double const f_s  = v[4];
    double const f_f  = 1.0 - f_s;
    double const R_fs = v[5] / f_f;
    double const R_sf = v[5] / f_s;
    double const B1   = v[6];

    // State vector is [x_f y_f z_f x_s z_s 1]
    AugMat       R        = Relax(f_f, R1_f, R2_f, f_s, R1_s, R2_s);
    AugMat       K        = Exchange(R_fs, R_sf);
    AugMat const RpK      = R + K;
    AugMat const S        = Spoil();
    AugMat const ramp     = S * (RpK * sequence.Tramp).exp();
    AugMat const spoilTRs = S * (RpK * sequence.TR * sequence.spoilers).exp();
    AugMat const Dprep    = (RpK * sequence.Dprep).exp();
    AugMat const Dseg     = (RpK * sequence.Dseg).exp();
    // QI_DBMAT(R)
    // QI_DBMAT(K)
    // QI_DBMAT(S)
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    std::vector<AugMat> prep_mats(sequence.preps());
    for (int ip = 0; ip < sequence.FA.rows(); ip++) {
        auto const rfa = RF(sequence.FA[ip], sequence.Trf, B1, R2_s, A_sl);
        A_mats[ip]     = ((RpK + rfa) * sequence.Trf).exp();
        R_mats[ip]     = (RpK * (sequence.TR - sequence.Trf)).exp();
        seg_mats[ip]   = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS) * spoilTRs;
        AugMat rfp     = RF(sequence.FAprep[ip], sequence.Tprep[ip], B1, R2_s, A_sl);
        prep_mats[ip]  = ((RpK + rfp) * sequence.Tprep[ip]).exp();
    }
    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int ip = 0; ip < sequence.preps(); ip++) {
        X = Dseg * ramp * seg_mats[ip] * ramp * Dprep * prep_mats[ip] * X;
    }
    AugVec const m_ss = SolveSteadyState(X);
    QI_DBVEC(m_ss);
    // Now loop through the segments and record the signal for each
    QI_ARRAY(DataType) sig(sequence.SPS * sequence.FAprep.size());
    AugVec m  = m_ss;
    int    ii = 0;
    for (int ip = 0; ip < sequence.preps(); ip++) {
        m = spoilTRs * ramp * Dprep * prep_mats[ip] * m;
        for (int is = 0; is < sequence.SPS; is++) {
            m         = A_mats[ip] * m;
            sig[ii++] = M0 * m[1];
            m         = S * R_mats[ip] * m;
        }
        m = Dseg * ramp * m;
    }
    if (basis.size()) {
        return basis * sig.matrix();
    } else {
        return sig;
    }
}

} // namespace QI
