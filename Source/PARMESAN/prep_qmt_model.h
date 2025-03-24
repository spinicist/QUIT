#pragma once

#include "GridInterp.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

namespace QI {
struct PrepQMTModel : Model<double, double, 4, 4, 1, 0, RealNoise<double>> {
    static int const NS = 1; // Number of parameters that need to be scaled
    PrepZTESequence const &sequence;
    RegularGrid const  &A_sl;
    double const T2_f = 0.07;
    double const T1_s = 0.35;
    double const T2_s = 14e-6;
    double const Rx = 14;
    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "f_s", "B1"};
    VaryingArray const                start{1.0, 2.0, 0.1, 1.0};
    VaryingArray const                lo{0.01, 0.5, 0.01, 0.5};
    VaryingArray const                hi{100.0, 5.0, 0.5, 1.5};

    std::array<std::string, NF> const fixed_names{"T2_f", "T1_s", "T2_s", "R_x"};
    FixedArray const                  fixed_defaults{T2_f, T1_s, T2_s, Rx};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    using PT     = ParameterType;
    using MT     = double;                   // Do the matrix maths at reduced precision for speed
    using AugMat = Eigen::Matrix<MT, 6, 6>; // Short for Augmented Matrix
    using AugVec = Eigen::Vector<MT, 6>;

    auto Relax(PT const M0_f,
               PT const R1_f,
               PT const R2_f,
               PT const M0_s,
               PT const R1_s,
               PT const R2_s) const -> AugMat {
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

    auto Exchange(PT const R_fs, PT const R_sf) const -> AugMat {
        AugMat K;
        K.setZero();
        K(2, 2) = -R_fs;
        K(2, 4) = R_sf;
        K(4, 2) = R_fs;
        K(4, 4) = -R_sf;
        return K;
    }

    auto RF(PT const flip, PT const pw, PT const B1, PT const R2_s) const -> AugMat {
        PT const B1x = B1 * flip / pw;
        PT const τ   = pw * R2_s;
        PT const ω   = B1 * flip / τ;
        AugMat  rf;
        rf.setZero();
        rf(1, 2) = B1x;
        rf(2, 1) = -B1x;
        rf(3, 3) = -R2_s * A_sl(ω, τ);
        rf(3, 4) = B1x;
        rf(4, 3) = -B1x;
        return rf;
    }

    auto Spoil() const -> AugMat {
        AugMat S;
        S.setZero();
        S(2, 2) = 1;
        S(4, 4) = 1;
        S(5, 5) = 1;
        return S;
    }

    auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType) {
        PT const M0   = v[0];
        PT const R1_f = 1. / v[1];
        PT const f_s  = v[2];
        PT const f_f  = 1.0 - f_s;
        PT const B1   = v[3];
        PT const R2_f = PT(1) / f[0];
        PT const R1_s = PT(1) / f[1];
        PT const R2_s = PT(1) / f[2];
        PT const Rx   = PT(f[3]);
        PT const M0_s = f_s * M0;
        PT const M0_f = M0 - M0_s;
        PT const R_fs = Rx * f_s;
        PT const R_sf = Rx * f_f;

        // State vector is [x_f y_f z_f x_s z_s 1]
        AugMat       R        = Relax(M0_f, R1_f, R2_f, M0_s, R1_s, R2_s);
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
            auto const rfa = RF(sequence.FA[ip], sequence.Trf, B1, R2_s);
            A_mats[ip]     = ((RpK + rfa) * sequence.Trf).exp();
            R_mats[ip]     = (RpK * (sequence.TR - sequence.Trf)).exp();
            seg_mats[ip]   = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS) * spoilTRs;
            AugMat rfp     = RF(sequence.FAprep[ip], sequence.Tprep, B1, R2_s);
            prep_mats[ip]  = ((RpK + rfp) * sequence.Tprep).exp();
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
                sig[ii++] = m[1];
                m         = S * R_mats[ip] * m;
            }
            m = Dseg * ramp * m;
        }
        if (sequence.basis.size()) {
            return sequence.basis * sig.matrix();
        } else {
            return sig;
        }
    }
};

    struct PrepQMTFullModel : Model<double, double, 7, 1, 1, 0, RealNoise<double>> {
        static int const NS = 1; // Number of parameters that need to be scaled
        PrepZTESequence const &sequence;
        RegularGrid const  &A_sl;
        // Pools are "free" and "semi-solid"
        std::array<std::string, NV> const varying_names{"M0", "T1_f", "T1_s", "T2_s", "f_s", "R_x", "B1"};
        VaryingArray const                start{1.0, 2.0, 1.0, 0.07, 0.35, 14e-6, 14.0};
        VaryingArray const                lo{0.01, 0.5, 0.5, 0.01, 0.05, 10e-6, 1.0};
        VaryingArray const                hi{100.0, 5.0, 1.5, 0.25, 0.5, 20e-6, 50.0};
    
        std::array<std::string, NF> const fixed_names{"T2_f"};
        FixedArray const                  fixed_defaults{0.1};

        int input_size(const int /* Unused */) const { return sequence.size(); }
    
        using PT     = ParameterType;
        using MT     = double;                   // Do the matrix maths at reduced precision for speed
        using AugMat = Eigen::Matrix<MT, 6, 6>; // Short for Augmented Matrix
        using AugVec = Eigen::Vector<MT, 6>;
    
        auto Relax(PT const M0_f,
                   PT const R1_f,
                   PT const R2_f,
                   PT const M0_s,
                   PT const R1_s,
                   PT const R2_s) const -> AugMat {
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
    
        auto Exchange(PT const R_fs, PT const R_sf) const -> AugMat {
            AugMat K;
            K.setZero();
            K(2, 2) = -R_fs;
            K(2, 4) = R_sf;
            K(4, 2) = R_fs;
            K(4, 4) = -R_sf;
            return K;
        }
    
        auto RF(PT const flip, PT const pw, PT const B1, PT const R2_s) const -> AugMat {
            PT const B1x = B1 * flip / pw;
            PT const τ   = pw * R2_s;
            PT const ω   = B1 * flip / τ;
            AugMat  rf;
            rf.setZero();
            rf(1, 2) = B1x;
            rf(2, 1) = -B1x;
            rf(3, 3) = -R2_s * A_sl(ω, τ);
            rf(3, 4) = B1x;
            rf(4, 3) = -B1x;
            return rf;
        }
    
        auto Spoil() const -> AugMat {
            AugMat S;
            S.setZero();
            S(2, 2) = 1;
            S(4, 4) = 1;
            S(5, 5) = 1;
            return S;
        }
    
        auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType) {
            // "M0", "T1_f", "T2_f", "T1_s", "T2_s", "f_s", "Rx", "B1"
            PT const M0   = v[0];
            PT const R1_f = 1. / v[1];
            PT const R2_f = PT(1) / f[0];
            PT const R1_s = PT(1) / v[2];
            PT const R2_s = PT(1) / v[3];
            PT const f_s  = v[4];
            PT const Rx   = v[5];
            PT const f_f  = 1.0 - f_s;
            PT const B1   = v[6];
            PT const M0_s = f_s * M0;
            PT const M0_f = M0 - M0_s;
            PT const R_fs = Rx * f_s;
            PT const R_sf = Rx * f_f;
    
            // State vector is [x_f y_f z_f x_s z_s 1]
            AugMat       R        = Relax(M0_f, R1_f, R2_f, M0_s, R1_s, R2_s);
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
                auto const rfa = RF(sequence.FA[ip], sequence.Trf, B1, R2_s);
                A_mats[ip]     = ((RpK + rfa) * sequence.Trf).exp();
                R_mats[ip]     = (RpK * (sequence.TR - sequence.Trf)).exp();
                seg_mats[ip]   = (S * R_mats[ip] * A_mats[ip]).pow(sequence.SPS) * spoilTRs;
                AugMat rfp     = RF(sequence.FAprep[ip], sequence.Tprep, B1, R2_s);
                prep_mats[ip]  = ((RpK + rfp) * sequence.Tprep).exp();
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
                    sig[ii++] = m[1];
                    m         = S * R_mats[ip] * m;
                }
                m = Dseg * ramp * m;
            }
            if (sequence.basis.size()) {
                return sequence.basis * sig.matrix();
            } else {
                return sig;
            }
        }
    };

} // namespace QI
