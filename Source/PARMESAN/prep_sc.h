#pragma once

#include "Macro.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

struct PrepModel : QI::Model<double, double, 5, 0, 1, 0, QI::RealNoise<double>> {
    static int const      NS = 1;
    PrepSequence const   &sequence;
    Eigen::MatrixXd const basis = Eigen::MatrixXd();
    static std::array<std::string, NV> const varying_names;
    static VaryingArray const start;
    static VaryingArray const lo;
    static VaryingArray const hi;

    std::array<std::string, NF> const fixed_names{};
    FixedArray const                  fixed_defaults{};

    auto input_size(const int /* Unused */) const -> int;
    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(DataType);

    auto dsdθ(VaryingArray const &v, FixedArray const &, int i) const -> QI_ARRAY(DataType);
};
