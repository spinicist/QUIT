#pragma once

#include "Model.h"
#include "transient_sequence.h"

struct MUPAModel : QI::Model<double, double, 3, 0> {
    static int const   NS = 1;
    RUFISSequence &    sequence;
    VaryingArray const start{30., 1., 0.1};
    VaryingArray const lo{1, 0.01, 0.01};
    VaryingArray const hi{150, 5.0, 5.0};

    std::array<std::string, NV> const varying_names{"M0", "T1", "T2"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
};

template <> struct QI::NoiseFromModelType<MUPAModel> : QI::RealNoise {};
