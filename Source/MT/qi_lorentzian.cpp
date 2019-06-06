/*
 *  qi_cestlorentz.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "JSON.h"
#include "MTSatSequence.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SequenceBase.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

template <int NP_> struct LorentzModel {
    using SequenceType  = QI::MTSatSequence;
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NP   = NP_;
    static constexpr int NVpP = 3; // Number of varying parameters per pool - A and FWHM
    static constexpr int NVg  = 0;
    static constexpr int NV   = NVg + NVpP * NP;
    static constexpr int ND   = 0;
    static constexpr int NF   = 0;

    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);

    SequenceType &              sequence;
    std::array<std::string, NV> varying_names;
    VaryingArray                bounds_lo;
    VaryingArray                bounds_hi;
    VaryingArray                start;
    std::array<bool, NP>        use_bandwidth;
    double                      Zref;
    bool                        additive;

    const std::array<std::string, NF> fixed_names{};
    const FixedArray                  fixed_defaults{};

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) &
                /* Unused */) const -> QI_ARRAY(typename Derived::Scalar) {
        using T       = typename Derived::Scalar;
        QI_ARRAY(T) S = QI_ARRAY(T)::Constant(sequence.sat_f0.rows(), T(Zref));
        QI_ARRAY(T) F;

        QI_ARRAY(T) const Z = QI_ARRAY(T)::Zero(sequence.sat_f0.rows());
        for (auto i = 0; i < NP; i++) {
            auto const indN = NVg + NVpP * i;
            T const &  df   = v[indN + 0];
            T const &  fwhm = v[indN + 1];
            T const &  A    = v[indN + 2];
            if (use_bandwidth[i]) {
                auto const x   = (sequence.sat_f0 - df - sequence.pulse.bandwidth / 2);
                auto const y   = (sequence.sat_f0 - df + sequence.pulse.bandwidth / 2);
                auto const xHx = (x > Z).select(x, Z);
                auto const yHy = (y < Z).select(y, Z);
                F              = xHx + yHy;
            } else {
                F = (sequence.sat_f0 - df - v[1]);
            }
            auto const L = A / (1.0 + (2.0 * F / fwhm).square());
            if (additive) {
                S += L;
            } else {
                S -= L;
            }
        }
        return S;
    }
};

// Declare here so available in Process function
args::ArgumentParser parser("Simple Lorentzian fitting.\nhttp://github.com/spinicist/QUIT");

args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
args::ValueFlag<int>
    pools(parser, "POOLS", "Number of Lorentzians to fit, default 1", {'p', "pools"}, 1);
QI_COMMON_ARGS;
args::Flag additive(
    parser, "ADDITIVE", "Use an additive model instead of subtractive", {'a', "add"}, false);
args::ValueFlag<double>
    Zref(parser, "Zref", "Reference value for Z-spectra, default 1.0", {'z', "zref"}, 1.0);

template <int N> void Process() {
    using LM   = LorentzModel<N>;
    using LFit = QI::NLLSFitFunction<LM>;

    QI::CheckPos(input_path);
    json input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::Log(verbose, "Reading sequence information");
    auto sequence = input.at("MTSat").get<QI::MTSatSequence>();

    auto const &pools_json = input.at("pools");
    if (pools_json.size() != N) {
        QI::Fail("Incorrect number of pools in JSON ({}, should be {})", pools_json.size(), N);
    }

    std::array<std::string, LM::NV> varying_names;
    typename LM::VaryingArray       low, high, start;
    std::array<bool, N>             use_bandwidth;
    for (size_t i = 0; i < N; i++) {
        auto const &pool = pools_json[i];
        auto const &name = pool["name"].get<std::string>();

        varying_names[LM::NVg + LM::NVpP * i + 0] = name + "_f0"s;
        varying_names[LM::NVg + LM::NVpP * i + 1] = name + "_fwhm"s;
        varying_names[LM::NVg + LM::NVpP * i + 2] = name + "_A"s;

        auto const &df_json   = pool["df0"];
        auto const &fwhm_json = pool["fwhm"];
        auto const &A_json    = pool["A"];
        if (df_json.size() != 3) {
            QI::Fail("Must specify start, low, high for df0 in {}", name);
        }
        if (fwhm_json.size() != 3) {
            QI::Fail("Must specify start, low, high for FWHM in {}", name);
        }
        if (A_json.size() != 3) {
            QI::Fail("Must specify start, low, high for A in {}", name);
        }

        start[LM::NVg + LM::NVpP * i + 0] = df_json[0].get<double>();
        start[LM::NVg + LM::NVpP * i + 1] = fwhm_json[0].get<double>();
        start[LM::NVg + LM::NVpP * i + 2] = A_json[0].get<double>();
        low[LM::NVg + LM::NVpP * i + 0]   = df_json[1].get<double>();
        low[LM::NVg + LM::NVpP * i + 1]   = fwhm_json[1].get<double>();
        low[LM::NVg + LM::NVpP * i + 2]   = A_json[1].get<double>();
        high[LM::NVg + LM::NVpP * i + 0]  = df_json[2].get<double>();
        high[LM::NVg + LM::NVpP * i + 1]  = fwhm_json[2].get<double>();
        high[LM::NVg + LM::NVpP * i + 2]  = A_json[2].get<double>();

        if (pool.find("use_bandwidth") != pool.end()) {
            use_bandwidth[i] = pool.at("use_bandwidth").get<bool>();
        } else {
            use_bandwidth[i] = false;
        }
    }

    QI::Log(verbose,
            "Fitting check:\nLow:   {}\nHigh:  {}\nStart:  {}\nBandwidth: {}",
            low.transpose(),
            high.transpose(),
            start.transpose(),
            fmt::join(use_bandwidth, ", "));

    LM model{sequence, varying_names, low, high, start, use_bandwidth, Zref.Get(), additive};

    if (simulate) {
        QI::SimulateModel<LM, false>(input, model, {}, {input_path.Get()}, verbose, simulate.Get());
    } else {
        LFit fit{model};
        auto fit_filter = QI::ModelFitFilter<LFit>::New(&fit, verbose, resids, subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, {}, mask.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "LTZ_");
        QI::Log(verbose, "Finished.");
    }
}

int main(int argc, char **argv) {
    Eigen::initParallel();

    QI::ParseArgs(parser, argc, argv, verbose, threads);
    switch (pools.Get()) {
    case 1:
        Process<1>();
        break;
    case 2:
        Process<2>();
        break;
    case 3:
        Process<3>();
        break;
    default:
        QI::Fail("Desired number of pools ({}) has not been implemented", pools.Get());
    }
    return EXIT_SUCCESS;
}