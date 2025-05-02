#pragma once

#include "GridInterp.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

namespace QI {
struct PrepQMTRx : Model<double, double, 4, 4, 1, 0, RealNoise<double>> {
    static int const      NS = 1; // Number of parameters that need to be scaled
    PrepSequence const   &sequence;
    InterpGrid const    &A_sl;
    Eigen::MatrixXd const basis = Eigen::MatrixXd();
    double const          T2_f  = 0.07;
    double const          T1_s  = 0.35;
    double const          T2_s  = 14e-6;
    double const          R_x   = 14;

    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "f_s", "B1"};
    VaryingArray const                start{1.0, 2.0, 0.1, 1.0};
    VaryingArray const                lo{0.01, 0.5, 0.01, 0.5};
    VaryingArray const                hi{100.0, 5.0, 0.5, 1.5};

    std::array<std::string, NF> const fixed_names{"T2_f", "T1_s", "T2_s", "R_x"};
    FixedArray const                  fixed_defaults{T2_f, T1_s, T2_s, R_x};

    int  input_size(const int /* Unused */) const;
    auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType);
};

struct PrepQMTRxFull : Model<double, double, 9, 0, 1, 0, RealNoise<double>> {
    static int const      NS = 1; // Number of parameters that need to be scaled
    PrepSequence const   &sequence;
    InterpGrid const    &A_sl;
    Eigen::MatrixXd const basis = Eigen::MatrixXd();
    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{
        "M0", "T1_f", "T2_f", "T1_s", "T2_s", "f_s", "R_x", "B1", "df0"};
    VaryingArray const start{1.0, 2.0, 0.07, 0.35, 14e-6, 0.1, 14.0, 1.0, 0.0};
    VaryingArray const lo{0.01, 0.5, 0.02, 0.05, 12e-6, 0.01, 1.0, 0.5, -150.0};
    VaryingArray const hi{100.0, 5.0, 2.5, 0.5, 20e-6, 0.5, 50.0, 1.5, 150.0};

    int  input_size(const int /* Unused */) const;
    auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType);
};

struct PrepQMTk : Model<double, double, 4, 4, 1, 0, RealNoise<double>> {
    static int const      NS = 1; // Number of parameters that need to be scaled
    PrepSequence const   &sequence;
    InterpGrid const    &A_sl;
    Eigen::MatrixXd const basis = Eigen::MatrixXd();
    double const          T2_f  = 0.07;
    double const          T1_s  = 0.35;
    double const          T2_s  = 14e-6;
    double const          k     = 1.4;

    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "f_s", "B1"};
    VaryingArray const                start{1.0, 2.0, 0.1, 1.0};
    VaryingArray const                lo{0.01, 0.5, 0.01, 0.5};
    VaryingArray const                hi{100.0, 5.0, 0.5, 1.5};

    std::array<std::string, NF> const fixed_names{"T2_f", "T1_s", "T2_s", "k"};
    FixedArray const                  fixed_defaults{T2_f, T1_s, T2_s, k};

    int  input_size(const int /* Unused */) const;
    auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType);
};

struct PrepQMTkFull : Model<double, double, 9, 1, 1, 0, RealNoise<double>> {
    static int const      NS = 1; // Number of parameters that need to be scaled
    PrepSequence const   &sequence;
    InterpGrid const    &A_sl;
    Eigen::MatrixXd const basis = Eigen::MatrixXd();
    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "T2_f", "T1_s", "T2_s", "f_s", "k", "B1", "df0"};
    VaryingArray const start{1.0, 2.0, 0.07, 0.35, 14e-6, 0.1, 1.4, 1.0, 0.0};
    VaryingArray const lo{0.01, 0.5, 0.02, 0.05, 12e-6, 0.01, 0.1, 0.5, -150.0};
    VaryingArray const hi{100.0, 5.0, 2.5, 0.5, 20e-6, 0.5, 10.0, 1.5, 150.0};

    std::array<std::string, NF> const fixed_names{"T2_f"};
    FixedArray const                  fixed_defaults{0.1};

    int  input_size(const int /* Unused */) const;
    auto signal(VaryingArray const &v, FixedArray const &f) const -> QI_ARRAY(DataType);
};

} // namespace QI
