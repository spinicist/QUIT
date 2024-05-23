#pragma once
#include "MTSequences.h"
#include "Model.h"

using namespace std::literals;

struct MTSatModel : QI::Model<double, double, 3, 1, 3, 0> {
    QI::MTSatSequence const &sequence;
    double const             C;
    bool const               smallangle;
    double const             delta_max;
    double const             r1_max;

    std::array<const std::string, NV> const varying_names{{"PD"s, "R1"s, "delta"s}};
    std::array<const std::string, 3> const  derived_names{};
    std::array<const std::string, NF> const fixed_names{{"B1"s}};

    FixedArray const fixed_defaults{1.0};

    int    input_size(int const /* Unused */) const { return 1; }
    int    output_size(int const /* Unused */) { return 1; }
    size_t num_outputs() const { return 3; }

    template <typename Derived>
    auto signals(const Eigen::ArrayBase<Derived> &v, const FixedArray &f) const
        -> std::array<QI_ARRAY(typename Derived::Scalar), NI> {
        auto const &S0    = v[0];
        auto const &R1    = v[1];
        auto const &dC    = v[2] / 100.;
        auto const &B1    = f[0];
        auto const &s     = sequence;
        auto        delta = dC * (1.0 - C * B1) / (1.0 - C);
        auto const  S_pd  = S0 * sin(B1 * s.al_pd) * (1. - exp(-R1 * s.TR_pd)) /
                          (1. - cos(B1 * s.al_pd) * exp(-R1 * s.TR_pd));
        auto const S_t1 = S0 * sin(B1 * s.al_t1) * (1. - exp(-R1 * s.TR_t1)) /
                          (1. - cos(B1 * s.al_t1) * exp(-R1 * s.TR_t1));
        auto const S_mt = S0 * sin(B1 * s.al_mt) * (1. - exp(-R1 * s.TR_mt)) /
                          (1. - (1. - delta) * cos(B1 * s.al_mt) * exp(-R1 * s.TR_mt));

        Eigen::ArrayXd a_pd(1), a_t1(1), a_mt(1);
        a_pd[0] = S_pd;
        a_t1[0] = S_t1;
        a_mt[0] = S_mt;

        return {a_pd, a_t1, a_mt};
    }

    QI_ARRAYN(ParameterType, NV)
    fit(QI_ARRAYN(DataType, NI) const &in, QI_ARRAYN(ParameterType, NF) const &fixed) const;
};
