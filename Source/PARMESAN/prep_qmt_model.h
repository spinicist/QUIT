#pragma once

#include "GridInterp.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

namespace QI {
struct PrepQMTModel : Model<double, double, 4, 4, 1, 0> {
    static int const NS = 1; // Number of parameters that need to be scaled
    PrepZTESequence &sequence;
    RegularGrid     &R2sl;
    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "f_s", "B1"};
    VaryingArray const                start{30.0, 2.0, 0.1, 1.0};
    VaryingArray const                lo{0.1, 0.5, 0.01, 0.5};
    VaryingArray const                hi{100.0, 5.0, 0.5, 1.5};

    std::array<std::string, NF> const fixed_names{"T2_f", "T1_s", "T2_s", "Rx"};
    FixedArray const                  fixed_defaults{0.07, 0.35, 14e-6, 14.0};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    template <typename Derived>
    auto signal(Eigen::ArrayBase<Derived> const &v, FixedArray const &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T      = typename Derived::Scalar;
        using ArrayT = Eigen::Array<T, -1, 1>;
        using AugMat = Eigen::Matrix<T, 6, 6>; // Short for Augmented Matrix
        using AugVec = Eigen::Vector<T, 6>;

        T const      M0   = v[0];
        T const      R1_f = 1. / v[1];
        T const      f_s  = v[2];
        T const      f_f  = 1.0 - f_s;
        T const      B1   = v[3];
        T const R2_f = T(1) / f[0];
        T const R1_s = T(1) / f[1];
        T const R2_s = T(1) / f[2];
    T const Rx   = T(f[3]);

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
        R << -T(R2_f), 0, 0, 0, 0, 0,             //
            0, -T(R2_f), 0, 0, 0, 0,              //
            0, 0, -R1_f, 0, 0, M0_f * R1_f,       //
            0, 0, 0, -T(R2_s), 0, 0,              //
            0, 0, 0, 0, -T(R1_s), M0_s * T(R1_s), //
            0, 0, 0, 0, 0, 0;
        AugMat K;
        K.setZero();
        K(2, 2) = -R_fs;
        K(2, 4) = R_sf;
        K(4, 2) = R_fs;
        K(4, 4) = -R_sf;
        AugMat const RpK = R + K;
        QI_DBMAT(R);
        QI_DBMAT(K);
        QI_DBMAT(RpK);
        using DT       = Eigen::DiagonalMatrix<T, 6, 6>;
        using DTV      = typename DT::DiagonalVectorType;
        AugMat const S = DT(DTV{T(0), T(0), T(1), T(0), T(1), T(1)}).toDenseMatrix();
        QI_DBMAT(S);
        // Setup readout segment matrices
        AugMat const ramp = S * (RpK * sequence.Tramp).exp();
        QI_DBMAT(ramp);
        std::vector<AugMat> A_mats(sequence.FA.rows());
        std::vector<AugMat> R_mats(sequence.FA.rows());
        std::vector<AugMat> seg_mats(sequence.FA.rows());
        std::vector<AugMat> prep_mats(sequence.preps());
        for (int ip = 0; ip < sequence.FA.rows(); ip++) {
            T const B1x  = B1 * sequence.FA[ip] / sequence.Trf;
            T const r2sl = T(0);
            AugMat  rf;
            rf.setZero();
            rf(1, 2) = B1x;
            rf(2, 1) = -B1x;
            rf(3, 3) = -R2_s;
            rf(3, 4) = B1x;
            rf(4, 3) = -B1x;
            A_mats[ip]   = ((RpK + rf) * sequence.Trf).exp();
            R_mats[ip]   = (RpK * (sequence.TR - sequence.Trf)).exp();
            seg_mats[ip] = (S * R_mats[ip] * A_mats[ip]).pow(T(sequence.SPS));

            T const B1p = B1 * sequence.FAprep[ip] / sequence.Tprep;
            AugMat  rfp;
            rfp.setZero();
            rf(1, 2) = B1p;
            rf(2, 1) = -B1p;
            rf(3, 3) = -R2_s;
            rf(3, 4) = B1p;
            rf(4, 3) = -B1p;
            QI_DBMAT(RpK);
            QI_DBMAT(rfp);
            QI_DB(sequence.Tprep);
            prep_mats[ip] = ((RpK + rfp) * sequence.Tprep).exp();
            QI_DBMAT(prep_mats[ip]);
        }
        // First calculate the system matrix
        AugMat X = AugMat::Identity();
        for (int ip = 0; ip < sequence.preps(); ip++) {
            X = ramp * seg_mats[ip] * ramp * prep_mats[ip] * X;
        }
        AugVec const m_ss = SolveSteadyState(X);

        // Now loop through the segments and record the signal for each
        ArrayT sig(sequence.SPS * sequence.FAprep.size());
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
};

template <> struct NoiseFromModelType<PrepQMTModel> : RealNoise {};
} // namespace QI
