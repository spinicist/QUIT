#pragma once

#include "ModelHelpers.h"
#include "mupa_sequence.h"

struct MUPAModel {
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 3; // Number of varying parameters
    static constexpr int ND = 0; // Number of derived parameters
    static constexpr int NF = 0; // Number of fixed parameters
    static constexpr int NI = 1; // Number of inputs

    using VaryingArray = QI_ARRAYN(ParameterType, NV); // Type for the varying parameter array
    using FixedArray   = QI_ARRAYN(ParameterType, NF); // Type for the fixed parameter array

    // Sequence paramter structs
    MUPASequence &sequence;

    // Fitting start point and bounds
    // The PD will be scaled by the fitting function to keep parameters roughly the same magnitude
    VaryingArray const start{30., 1., 0.1};
    VaryingArray const lo{1, 0.01, 0.01};
    VaryingArray const hi{150, 5.0, 5.0};

    std::array<std::string, NV> const varying_names{"PD", "T1", "T2"};
    std::array<std::string, NF> const fixed_names{};
    // If fixed parameters not supplied, use these default values
    FixedArray const fixed_defaults{};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
};

template <> struct QI::NoiseFromModelType<MUPAModel> : QI::RealNoise {};

Eigen::Matrix<double, 4, 4> RF_1c(double const &B1x, double const &B1y);