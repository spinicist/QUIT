#pragma once

#include "Lineshape.h"
#include "Model.h"
#include "ss_sequence.h"

struct SS_MT_Model : QI::Model<double, double, 8, 0, 1, 1> {
    static int const    NS = 2;
    SSSequence &        sequence;
    QI::InterpLineshape lineshape;
    VaryingArray const  start{30.0, 3.0, 1.0, 0.1, 12e-6, 30., 0., 1.0};
    VaryingArray const  lo{0.1, 5e-6, 0.5, 0.005, 5e-6, 1., -250., 0.5};
    VaryingArray const  hi{60.0, 30.0, 5.0, 5.0, 25e-6, 100., 250., 1.5};

    std::array<std::string, NV> const varying_names{
        "M0_f", "M0_b", "T1_f", "T2_f", "T2_b", "k", "f0", "B1"};
    std::array<std::string, ND> const derived_names{"f_b"};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(double);
    void derived(const VaryingArray &varying, const FixedArray &, DerivedArray &derived) const;
};

template <> struct QI::NoiseFromModelType<SS_MT_Model> : QI::RealNoise {};
