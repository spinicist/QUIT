#pragma once

#include "GridInterp.h"
#include "Model.h"
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

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
    void derived(const VaryingArray &varying, const FixedArray &, DerivedArray &derived) const;
};

template <> struct NoiseFromModelType<PrepQMTModel> : RealNoise {};
} // namespace QI
