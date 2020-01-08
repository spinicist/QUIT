#pragma once

#include "Lineshape.h"
#include "Model.h"
#include "mupa_sequence.h"

struct MUPAMTB1Model : QI::Model<double, double, 4, 1, 1, 1> {
    MUPASequence &      sequence;
    QI::InterpLineshape lineshape;
    VaryingArray const  start{30.0, 1.0, 0.1, 3.0};
    VaryingArray const  lo{0.1, 0.5, 0.005, 0.1};
    VaryingArray const  hi{60.0, 5.0, 5.0, 20.0};

    VaryingNames const varying_names{"M0_f", "T1_f", "T2_f", "M0_b"};
    FixedNames const   fixed_names{"B1"};
    FixedArray const   fixed_defaults{1.0};
    DerivedNames const derived_names{"f_b"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
    void derived(const VaryingArray &varying, const FixedArray &, DerivedArray &derived) const;
};

template <> struct QI::NoiseFromModelType<MUPAMTB1Model> : QI::RealNoise {};
