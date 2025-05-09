#pragma once

#include "Macro.h"
#include "Model.h"
#include "parmesan.hpp"
#include "prep_sequence.h"

struct PrepModel : QI::Model<double, double, 5, 0, 1, 0, QI::RealNoise<double>> {
    static int const                  NS = 1;
    PrepSequence const               &sequence;
    Eigen::MatrixXd const             basis = Eigen::MatrixXd();
    std::array<std::string, NV> const varying_names{"M0", "T1", "T2", "B1", "df"};
    VaryingArray const                start{1., 1., 0.1, 1.0, 0.0};
    VaryingArray const                lo{0.01, 0.01, 0.01, 0.5, -250.0};
    VaryingArray const                hi{1000., 5.0, 3.5, 1.5, 260.0};

    std::array<std::string, NF> const fixed_names{};
    FixedArray const                  fixed_defaults{};

    auto                              input_size(const int /* Unused */) const -> int;
    auto signal(VaryingArray const &v, FixedArray const &) const -> QI_ARRAY(DataType);

    auto dsdθ(VaryingArray const &v, FixedArray const &, int i) const -> QI_ARRAY(DataType);
};
