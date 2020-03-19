#pragma once

#include "Lineshape.h"
#include "Model.h"
#include "mt_sequence.h"

struct MTModel : QI::Model<double, double, 7, 0, 1, 1> {
    MTSequence &        sequence;
    QI::InterpLineshape lineshape;
    VaryingArray const  start{30.0, 3.0, 1.0, 0.1, 12e-6, 5., 1.0};
    VaryingArray const  lo{0.1, 5e-6, 0.5, 0.005, 5e-6, 1., 0.5};
    VaryingArray const  hi{60.0, 30.0, 5.0, 5.0, 25e-6, 10., 1.5};

    std::array<std::string, NV> const varying_names{
        "M0_f", "M0_b", "T1_f", "T2_f", "T2_b", "k", "B1"};
    std::array<std::string, ND> const derived_names{"f_b"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
    void derived(const VaryingArray &varying, const FixedArray &, DerivedArray &derived) const;
};
