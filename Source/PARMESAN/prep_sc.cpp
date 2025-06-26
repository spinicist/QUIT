// #define QI_DEBUG_BUILD 1
#include "prep_sc.h"
#include "Log.h"
namespace {

using AugMat = Eigen::Matrix<double, 4, 4>;
using AugVec = Eigen::Vector<double, 4>;

auto Relax(double const M0, double const R1, double const R2, double const df0) -> AugMat {
    AugMat R;
    R.setZero();
    double const dω0 = df0 * 2.0 * M_PI;
    R(0, 0)          = -R2;
    R(0, 1)          = dω0;
    R(1, 0)          = -dω0;
    R(1, 1)          = -R2;
    R(2, 2)          = -R1;
    R(2, 3)          = R1 * M0;
    return R;
}

auto RF(double const ω0, double const ω1x) -> AugMat {
    AugMat rf;
    rf.setZero();
    rf(0, 1) = ω0;
    rf(1, 0) = -ω0;
    rf(1, 2) = ω1x;
    rf(2, 1) = -ω1x;
    return rf;
}

auto BlockPulse(double const flip, double const pw, double const f0, AugMat const R) -> AugMat {
    double const ω1 = flip / pw;
    double const ω0 = f0 * 2.0 * M_PI;
    return ((RF(ω0, ω1) + R) * pw).exp();
}

auto GaussPulse(double const flip, double const pw, double const f0, AugMat const R) -> AugMat {
    double const sigma = 0.3;
    double const ω0    = f0 * 2.0 * M_PI;
    AugMat       gp    = AugMat::Identity();
    int const    N     = 64;
    double const tau   = pw / N;

    // Get the pulse area
    double a = 0.0;
    for (int ii = 0; ii < N; ii++) {
        double const t = tau * ii;
        a += std::exp(-0.5 * std::pow((t - (pw / 2.0)) / (sigma * (pw / 2.0)), 2)) * tau;
    }

    for (int ii = 0; ii < N; ii++) {
        double const t = tau * ii;
        double const ω1 =
            std::exp(-0.5 * std::pow((t - (pw / 2.0)) / (sigma * (pw / 2.0)), 2)) * flip / a;
        QI_DB(ω1);
        gp = ((RF(ω0, ω1) + R) * tau).exp() * gp;
    }
    QI_DB(tau);
    QI_DB(flip);
    QI_DB(pw);
    QI_DB(a);
    QI_DB(f0);
    QI_DBMAT(R);
    QI_DBMAT(gp);
    return gp;
}

auto Spoil() -> AugMat {
    AugMat S;
    S.setZero();
    S(2, 2) = 1;
    S(3, 3) = 1;
    return S;
}

} // namespace

std::array<std::string, PrepModel::NV> const PrepModel::varying_names{"M0", "T1", "T2", "B1", "df"};
PrepModel::VaryingArray const                PrepModel::start{1., 1., 0.1, 1.0, 0.0};
PrepModel::VaryingArray const                PrepModel::lo{0.01, 0.01, 0.01, 0.5, -250.0};
PrepModel::VaryingArray const                PrepModel::hi{1000., 5.0, 3.5, 1.5, 260.0};
std::array<std::string, PrepModel::NF> const PrepModel::fixed_names{};
PrepModel::FixedArray const                  PrepModel::fixed_defaults{};

auto PrepModel::input_size(const int /* Unused */) const -> int {
    if (basis.size()) {
        return basis.rows();
    } else {
        return sequence.size();
    }
}

auto PrepModel::signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(DataType) {

    double const M0  = v[0];
    double const R1  = 1. / v[1];
    double const R2  = 1. / v[2];
    double const B1  = v[3];
    double const df0 = v[4];

    QI_DBVEC(v);
    QI_DB(M0);
    QI_DB(R1);
    QI_DB(R2);
    QI_DB(B1);
    QI_DB(df0);
    QI_DB(sequence.Trf)
    QI_DB(sequence.Tprep)
    QI_DBVEC(sequence.FA * 180 / M_PI)
    QI_DBVEC(sequence.FAprep * 180 / M_PI)

    AugMat const R        = Relax(1.0, R1, R2, df0);
    AugMat const S        = Spoil();
    AugMat const ramp     = S * (R * sequence.Tramp).exp();
    AugMat const spoilTRs = S * (R * sequence.TR * sequence.spoilers).exp();

    QI_DBMAT(R);
    QI_DBMAT(S);
    QI_DBMAT(ramp);
    QI_DBMAT(spoilTRs);
    // Setup pulse matrices
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    std::vector<AugMat> prep_mats(sequence.preps());
    std::vector<AugMat> pre_mats(sequence.preps());
    std::vector<AugMat> post_mats(sequence.preps());
    for (int ip = 0; ip < sequence.preps(); ip++) {
        A_mats[ip]    = BlockPulse(sequence.FA[ip] * B1, sequence.Trf, 0, R);
        R_mats[ip]    = (R * (sequence.TR - sequence.Trf)).exp();
        pre_mats[ip]  = (R * sequence.Tpreseg[ip]).exp();
        post_mats[ip] = (R * sequence.Tpostseg[ip]).exp();
        seg_mats[ip]  = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS) * spoilTRs;
        QI_DBMAT(A_mats[ip]);
        QI_DBMAT(R_mats[ip]);
        QI_DBMAT(seg_mats[ip]);
        prep_mats[ip] =
            gauss ?
                GaussPulse(sequence.FAprep[ip] * B1, sequence.Tprep[ip], sequence.fprep[ip], R) :
                BlockPulse(sequence.FAprep[ip] * B1, sequence.Tprep[ip], sequence.fprep[ip], R);
        QI_DBMAT(prep_mats[ip]);
    }
    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int ip = 0; ip < sequence.preps(); ip++) {
        QI_DBMAT(X);
        X = post_mats[ip] * ramp * seg_mats[ip] * ramp * pre_mats[ip] * prep_mats[ip] * X;
    }
    AugVec const m_ss = SolveSteadyState(X);
    QI_DBMAT(X);
    QI_DBVEC(m_ss);
    // Now loop through the segments and record the signal for each
    QI_ARRAY(DataType) sig(sequence.SPS * sequence.FAprep.size());

    AugVec m  = m_ss;
    int    ii = 0;
    for (int ip = 0; ip < sequence.preps(); ip++) {
        QI_DBVEC(prep_mats[ip] * m);
        QI_DBVEC(ramp * prep_mats[ip] * m);
        m = spoilTRs * ramp * pre_mats[ip] * prep_mats[ip] * m;
        QI_DBVEC(m);
        for (int is = 0; is < sequence.SPS; is++) {
            m         = A_mats[ip] * m;
            sig[ii++] = M0 * m[1];
            m         = S * R_mats[ip] * m;
        }
        m = post_mats[ip] * ramp * m;
    }
    QI_DBVEC(sig);
    // exit(EXIT_FAILURE);
    if (basis.size()) {
        if (basis.cols() != sig.size()) {
            QI::Fail("Basis had {} columns, should be {}", basis.cols(), sig.size());
        }
        return basis * sig.matrix();
    } else {
        return sig;
    }
}

auto PrepModel::dsdθ(VaryingArray const &v, FixedArray const &f, int i) const
    -> QI_ARRAY(DataType) {
    // Central differences
    double const h = std::max(
        std::abs(v(i)) * 1e-8,
        1e-8); // From second answer on
               // https://math.stackexchange.com/questions/815113/is-there-a-general-formula-for-estimating-the-step-size-h-in-numerical-different
    auto vph = v, vmh = v;
    vph(i) += h;
    vmh(i) -= h;
    QI_ARRAY(DataType) const d = (signal(vph, f) - signal(vmh, f)) / (2 * h);
    QI_DBVEC(d);
    return d;
}

std::array<std::string, PrepModel2::NV> const PrepModel2::varying_names{"M0", "T1", "T2", "B1"};
PrepModel2::VaryingArray const                PrepModel2::start{1., 1., 0.1, 1.0};
PrepModel2::VaryingArray const                PrepModel2::lo{0.01, 0.01, 0.01, 0.5};
PrepModel2::VaryingArray const                PrepModel2::hi{1000., 5.0, 3.5, 1.5};
std::array<std::string, PrepModel2::NF> const PrepModel2::fixed_names{"df"};
PrepModel2::FixedArray const                  PrepModel2::fixed_defaults{0.0};

auto PrepModel2::input_size(const int /* Unused */) const -> int {
    if (basis.size()) {
        return basis.rows();
    } else {
        return sequence.size();
    }
}

auto PrepModel2::signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType) {

    double const M0  = v[0];
    double const R1  = 1. / v[1];
    double const R2  = 1. / v[2];
    double const B1  = v[3];
    double const df0 = f[0];

    QI_DBVEC(v);
    QI_DB(M0);
    QI_DB(R1);
    QI_DB(R2);
    QI_DB(B1);
    QI_DB(df0);
    QI_DB(sequence.Trf)
    QI_DB(sequence.Tprep)
    QI_DBVEC(sequence.FA * 180 / M_PI)
    QI_DBVEC(sequence.FAprep * 180 / M_PI)

    AugMat const R        = Relax(1.0, R1, R2, df0);
    AugMat const S        = Spoil();
    AugMat const ramp     = S * (R * sequence.Tramp).exp();
    AugMat const spoilTRs = S * (R * sequence.TR * sequence.spoilers).exp();

    QI_DBMAT(R);
    QI_DBMAT(S);
    QI_DBMAT(ramp);
    QI_DBMAT(spoilTRs);
    // Setup pulse matrices
    std::vector<AugMat> A_mats(sequence.FA.rows());
    std::vector<AugMat> R_mats(sequence.FA.rows());
    std::vector<AugMat> seg_mats(sequence.FA.rows());
    std::vector<AugMat> prep_mats(sequence.preps());
    std::vector<AugMat> pre_mats(sequence.preps());
    std::vector<AugMat> post_mats(sequence.preps());
    for (int ip = 0; ip < sequence.preps(); ip++) {
        A_mats[ip]    = BlockPulse(sequence.FA[ip] * B1, sequence.Trf, 0, R);
        R_mats[ip]    = (R * (sequence.TR - sequence.Trf)).exp();
        pre_mats[ip]  = (R * sequence.Tpreseg[ip]).exp();
        post_mats[ip] = (R * sequence.Tpostseg[ip]).exp();
        seg_mats[ip]  = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS) * spoilTRs;
        QI_DBMAT(A_mats[ip]);
        QI_DBMAT(R_mats[ip]);
        QI_DBMAT(seg_mats[ip]);
        prep_mats[ip] =
            GaussPulse(sequence.FAprep[ip] * B1, sequence.Tprep[ip], sequence.fprep[ip], R);
        QI_DBMAT(prep_mats[ip]);
    }
    // First calculate the system matrix
    AugMat X = AugMat::Identity();
    for (int ip = 0; ip < sequence.preps(); ip++) {
        QI_DBMAT(X);
        X = post_mats[ip] * ramp * seg_mats[ip] * ramp * pre_mats[ip] * prep_mats[ip] * X;
    }
    AugVec const m_ss = SolveSteadyState(X);
    QI_DBMAT(X);
    QI_DBVEC(m_ss);
    // Now loop through the segments and record the signal for each
    QI_ARRAY(DataType) sig(sequence.SPS * sequence.FAprep.size());

    AugVec m  = m_ss;
    int    ii = 0;
    for (int ip = 0; ip < sequence.preps(); ip++) {
        QI_DBVEC(prep_mats[ip] * m);
        QI_DBVEC(ramp * prep_mats[ip] * m);
        m = spoilTRs * ramp * pre_mats[ip] * prep_mats[ip] * m;
        QI_DBVEC(m);
        for (int is = 0; is < sequence.SPS; is++) {
            m         = A_mats[ip] * m;
            sig[ii++] = M0 * m[1];
            m         = S * R_mats[ip] * m;
        }
        m = post_mats[ip] * ramp * m;
    }
    QI_DBVEC(sig);
    // exit(EXIT_FAILURE);
    if (basis.size()) {
        if (basis.cols() != sig.size()) {
            QI::Fail("Basis had {} columns, should be {}", basis.cols(), sig.size());
        }
        return basis * sig.matrix();
    } else {
        return sig;
    }
}

auto PrepModel2::dsdθ(VaryingArray const &v, FixedArray const &f, int i) const
    -> QI_ARRAY(DataType) {
    // Central differences
    double const h = std::max(
        std::abs(v(i)) * 1e-8,
        1e-8); // From second answer on
               // https://math.stackexchange.com/questions/815113/is-there-a-general-formula-for-estimating-the-step-size-h-in-numerical-different
    auto vph = v, vmh = v;
    vph(i) += h;
    vmh(i) -= h;
    QI_ARRAY(DataType) const d = (signal(vph, f) - signal(vmh, f)) / (2 * h);
    QI_DBVEC(d);
    return d;
}
