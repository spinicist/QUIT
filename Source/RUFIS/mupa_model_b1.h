#pragma once

#include "Model.h"
#include "mupa_sequence.h"

struct MUPAB1Model : QI::Model<double, double, 4, 0> {
    static int const   NS = 1;
    MUPASequence &     sequence;
    VaryingArray const start{30., 1., 0.1, 1.0};
    VaryingArray const lo{1, 0.01, 0.01, 0.5};
    VaryingArray const hi{150, 5.0, 5.0, 1.5};

    std::array<std::string, NV> const varying_names{"M0", "T1", "T2", "B1"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
};

template <> struct QI::NoiseFromModelType<MUPAB1Model> : QI::RealNoise {};
