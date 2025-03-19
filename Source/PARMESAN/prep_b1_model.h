#pragma once

#include "Model.h"
#include "prep_sequence.h"

struct PrepB1Model : QI::Model<double, double, 3, 0> {
    static int const   NS = 1;
    PrepZTESequence   &sequence;
    std::array<std::string, NV> const varying_names{"M0", "T1", "B1"};
    VaryingArray const start{30., 1., 1.0};
    VaryingArray const lo{1, 0.01, 0.5};
    VaryingArray const hi{150, 5.0, 1.5};
    int input_size(const int /* Unused */) const { return sequence.size(); }
    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
};
template <> struct QI::NoiseFromModelType<PrepB1Model> : QI::RealNoise {};
