#pragma once

#include "Model.h"
#include "mupa_sequence.h"

struct MUPAMTModel : QI::Model<double, double, 5, 0, 1, 1> {
    static int const   NS = 2;
    MUPASequence &     sequence;
    VaryingArray const start{30.0, 3.0, 1.0, 0.1, 1.0};
    VaryingArray const lo{0.1, 0.1, 0.5, 0.005, 0.5};
    VaryingArray const hi{60.0, 60.0, 5.0, 5.0, 1.5};

    std::array<std::string, NV> const varying_names{"M0_f", "M0_b", "T1_f", "T2_f", "B1"};
    std::array<std::string, ND> const derived_names{"f_b"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
    void derived(const VaryingArray &varying, const FixedArray &, DerivedArray &derived) const;
};

template <> struct QI::NoiseFromModelType<MUPAMTModel> : QI::RealNoise {};
