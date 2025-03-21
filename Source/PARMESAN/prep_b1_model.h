#pragma once

#include "Macro.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

struct PrepB1Model : QI::Model<double, double, 3, 0> {
    static int const                  NS = 1;
    PrepZTESequence                  &sequence;
    std::array<std::string, NV> const varying_names{"M0", "T1", "B1"};
    VaryingArray const                start{1., 1., 1.0};
    VaryingArray const                lo{0.1, 0.01, 0.5};
    VaryingArray const                hi{100., 5.0, 1.5};
    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double) {
        using T      = double;
        using AugMat = Eigen::Matrix<T, 3, 3>;
        using AugVec = Eigen::Vector<T, 3>;

        T const      M0 = v[0];
        T const      R1 = 1. / v[1];
        double const R2 = 10.;
        T const      B1 = v[2];

        QI_DBVEC(v);
        // QI_DB(sequence.Trf)
        // QI_DB(sequence.Tprep)
        // QI_DBVEC(sequence.FA * 180 / M_PI)
        // QI_DBVEC(sequence.FAprep * 180 / M_PI)

        AugMat R;
        R << -R2, 0, 0,      //
            0, -R1, R1 * M0, //
            0, 0, 0;
        // QI_DBMAT(R);
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
            prep_mats[ip] = ((R + rfp) * sequence.Tprep).exp();
            // QI_DBMAT((rfp * sequence.Tprep).exp());
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
            m = ramp * prep_mats[ip] * m;
            for (int is = 0; is < sequence.SPS; is++) {
                m         = A_mats[ip] * m;
                sig[ii++] = m[1];
                m         = S * R_mats[ip] * m;
            }
            m = ramp * m;
        }
        // QI_DBVEC(sig);
        if (sequence.basis.size()) {
            return sequence.basis * sig.matrix();
        } else {
            return sig;
        }
    }
};
template <> struct QI::NoiseFromModelType<PrepB1Model> : QI::RealNoise {};
