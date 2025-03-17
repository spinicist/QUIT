#pragma once

#include "Model.h"
#include "transient_sequence.h"

struct QMTModel : QI::Model<double, double, 4, 4, 1, 1> {
    static int const   NS = 1;
    PrepSequence &    sequence;
    VaryingArray const start{30.0, 2.0, 0.1, 1.0};
    VaryingArray const lo{0.1, 0.5, 0.01, 0.5};
    VaryingArray const hi{100.0, 5.0, 0.5, 1.5};

    // Pools are "free" and "semi-solid"
    std::array<std::string, NV> const varying_names{"M0", "T1_f", "f_s", "B1"};
    std::array<std::string, ND> const fixed_names{"T2_f", "T1_s", "T2_s", "k"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
    void derived(const VaryingArray &varying, const FixedArray &, DerivedArray &derived) const;
};

template <> struct QI::NoiseFromModelType<MUPAMTModel> : QI::RealNoise {};
