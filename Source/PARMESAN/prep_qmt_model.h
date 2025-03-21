#pragma once

#include "GridInterp.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

namespace QI {
struct PrepQMTModel : Model<double, double, 4, 4, 1, 0> {
    static int const NS = 1; // Number of parameters that need to be scaled
    PrepZTESequence &sequence;
    RegularGrid     &A_sl;
    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "f_s", "B1"};
    VaryingArray const                start{30.0, 2.0, 0.1, 1.0};
    VaryingArray const                lo{0.1, 0.5, 0.01, 0.5};
    VaryingArray const                hi{100.0, 5.0, 0.5, 1.5};

    std::array<std::string, NF> const fixed_names{"T2_f", "T1_s", "T2_s", "Rx"};
    FixedArray const                  fixed_defaults{0.07, 0.35, 14e-6, 14.0};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    using T      = double;
    using AugMat = Eigen::Matrix<T, 6, 6>; // Short for Augmented Matrix
    using AugVec = Eigen::Vector<T, 6>;

    auto
    Relax(T const M0_f, T const R1_f, T const R2_f, T const M0_s, T const R1_s, T const R2_s) const
        -> AugMat {
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

    auto Exchange(T const R_fs, T const R_sf) const -> AugMat {
        AugMat K;
        K.setZero();
        K(2, 2) = -R_fs;
        K(2, 4) = R_sf;
        K(4, 2) = R_fs;
        K(4, 4) = -R_sf;
        return K;
    }

    auto RF(T const flip, T const pw, T const B1, T const R2_s) const -> AugMat {
        T const B1x = B1 * flip / pw;
        T const τ   = pw * R2_s;
        T const ω   = B1 * flip / τ;
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

    auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(T) {
        T const M0   = v[0];
        T const R1_f = 1. / v[1];
        T const f_s  = v[2];
        T const f_f  = 1.0 - f_s;
        T const B1   = v[3];
        T const R2_f = T(1) / f[0];
        T const R1_s = T(1) / f[1];
        T const R2_s = T(1) / f[2];
        T const Rx   = T(f[3]);
        T const M0_s = f_s * M0;
        T const M0_f = M0 - M0_s;
        T const R_fs = Rx * f_s;
        T const R_sf = Rx * f_f;

        // State vector is [x_f y_f z_f x_s z_s 1]
        AugMat              R    = Relax(M0_f, R1_f, R2_f, M0_s, R1_s, R2_s);
        AugMat              K    = Exchange(R_fs, R_sf);
        AugMat const        RpK  = R + K;
        AugMat const        S    = Spoil();
        AugMat const        ramp = S * (RpK * sequence.Tramp).exp();
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
            seg_mats[ip]   = (S * R_mats[ip] * A_mats[ip]).pow(T(sequence.SPS));
            AugMat rfp     = RF(sequence.FAprep[ip], sequence.Tprep, B1, R2_s);
            prep_mats[ip]  = ((RpK + rfp) * sequence.Tprep).exp();
        }
        // First calculate the system matrix
        AugMat X = AugMat::Identity();
        for (int ip = 0; ip < sequence.preps(); ip++) {
            X = ramp * seg_mats[ip] * ramp * prep_mats[ip] * X;
        }
        AugVec const m_ss = SolveSteadyState(X);
        QI_DBVEC(m_ss);
        // Now loop through the segments and record the signal for each
        QI_ARRAY(T) sig(sequence.SPS * sequence.FAprep.size());
        AugVec m  = m_ss;
        int    ii = 0;
        for (int ip = 0; ip < sequence.preps(); ip++) {
            m = ramp * prep_mats[ip] * m;
            for (int is = 0; is < sequence.SPS; is++) {
                m         = A_mats[ip] * m;
                sig[ii++] = m[1];
                m         = S * R_mats[ip] * m;
            }
            m = ramp * m;
        }
        if (sequence.basis.size()) {
            return sequence.basis * sig.matrix();
        } else {
            return sig;
        }
    }
};

template <> struct NoiseFromModelType<PrepQMTModel> : RealNoise {};
} // namespace QI
