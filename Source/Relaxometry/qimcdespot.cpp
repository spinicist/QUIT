/*
 *  qmcdespot.cpp
 *
 *  Created by Tobias Wood on 03/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>
#include <array>

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "RegionContraction.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "ThreePoolModel.h"
#include "TwoPoolModel.h"
#include "Util.h"

template <typename Model> struct MCDSRCFunctor {
    const Eigen::ArrayXd data, weights;
    const QI_ARRAYN(double, Model::NF) fixed;
    const Model &model;

    MCDSRCFunctor(const Model &m,
                  const QI_ARRAYN(double, Model::NF) & f,
                  const Eigen::ArrayXd &d,
                  const Eigen::ArrayXd &w) :
        data(d),
        weights(w), fixed(f), model(m) {}

    int inputs() const { return Model::NV; }
    int values() const { return model.spgr.size() + model.ssfp.size(); }

    bool constraint(const QI_ARRAYN(double, Model::NV) & varying) const {
        return model.valid(varying);
    }

    Eigen::ArrayXd residuals(const QI_ARRAYN(double, Model::NV) & varying) const {
        return data - model.signal(varying, fixed);
    }

    double operator()(const QI_ARRAYN(double, Model::NV) & varying) const {
        return (residuals(varying) * weights).square().sum();
    }
};

template <typename Model> struct SRCFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using RMSErrorType        = double;
    using FlagType            = int;
    using ModelType           = Model;
    Model &model;

    int input_size(const int &i) const {
        if (i == 0) {
            return model.spgr.size();
        } else if (i == 1) {
            return model.ssfp.size();
        } else {
            QI::Fail("Incorrect input number {}", i);
        }
    }
    int n_outputs() const { return Model::NV; }

    int    max_iterations = 5;
    size_t src_samples = 5000, src_retain = 50;
    bool   src_gauss = true;

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          typename Model::FixedArray const & fixed,
                          typename Model::VaryingArray &     v,
                          typename Model::CovarArray * /*Unused */,
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const {
        Eigen::ArrayXd data(model.ssfp.size() + model.spgr.size());
        int            dataIndex = 0;
        for (size_t i = 0; i < inputs.size(); i++) {
            if (model.scale_to_mean) {
                data.segment(dataIndex, inputs[i].rows()) = inputs[i] / inputs[i].mean();
            } else {
                data.segment(dataIndex, inputs[i].rows()) = inputs[i];
            }
            dataIndex += inputs[i].rows();
        }
        QI_ARRAYN(double, Model::NV) thresh = QI_ARRAYN(double, Model::NV)::Constant(0.05);
        const double & f0                   = fixed[0];
        Eigen::ArrayXd weights(model.spgr.size() + model.ssfp.size());
        weights.head(model.spgr.size()) = 1;
        weights.tail(model.ssfp.size()) = model.ssfp.weights(f0);
        using Functor                   = MCDSRCFunctor<Model>;
        Functor                        func(model, fixed, data, weights);
        QI::RegionContraction<Functor> rc(func,
                                          model.bounds_lo,
                                          model.bounds_hi,
                                          thresh,
                                          src_samples,
                                          src_retain,
                                          max_iterations,
                                          0.02,
                                          src_gauss,
                                          false);
        if (!rc.optimise(v)) {
            return {false, "Region contraction failed"};
        }
        auto r   = func.residuals(v);
        residual = sqrt(r.square().sum() / r.rows());
        if (residuals.size() > 0) {
            residuals[0] = r.head(model.spgr.size());
            residuals[1] = r.tail(model.ssfp.size());
        }
        iterations = rc.contractions();
        return {true, ""};
    }
};

//******************************************************************************
// Main
//******************************************************************************
int mcdespot_main(args::Subparser &parser) {
    args::Positional<std::string> spgr_path(parser, "SPGR FILE", "Input SPGR file");
    args::Positional<std::string> ssfp_path(parser, "SSFP FILE", "Input SSFP file");
    QI_COMMON_ARGS;
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (Hertz)", {'f', "f0"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<int>         modelarg(
        parser, "MODEL", "Select model to fit - 2/3, default 3", {'M', "model"}, 3);
    args::Flag scale(parser, "SCALE", "Normalize signals to mean (a good idea)", {'S', "scale"});
    args::Flag use_src(
        parser, "SRC", "Use flat prior (stochastic region contraction), not gaussian", {"SRC"});
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations, default 4", {'i', "its"}, 4);
    args::Flag           bounds(parser, "BOUNDS", "Specify bounds in input", {"bounds"});
    parser.Parse();
    QI::CheckPos(spgr_path);
    QI::CheckPos(ssfp_path);

    QI::Log(verbose, "Reading sequences");
    auto input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto spgr  = input.at("SPGR").get<QI::SPGRSequence>();
    auto ssfp  = input.at("SSFP").get<QI::SSFPSequence>();

    auto process = [&](auto model, const std::string &model_name) {
        if (simulate) {
            QI::SimulateModel<decltype(model), true>(input,
                                                     model,
                                                     {f0.Get(), B1.Get()},
                                                     {spgr_path.Get(), ssfp_path.Get()},
                                                     mask.Get(),
                                                     verbose,
                                                     simulate.Get(),
                                                     subregion.Get());
        } else {
            using FitType = SRCFit<decltype(model)>;
            FitType src{model};
            src.src_gauss = !use_src;
            if (bounds) {
                src.model.bounds_lo = QI::ArrayFromJSON<double>(input, "lower_bounds");
                src.model.bounds_hi = QI::ArrayFromJSON<double>(input, "upper_bounds");
            }
            QI::Log(verbose, "Low bounds: {}", src.model.bounds_lo.transpose());
            QI::Log(verbose, "High bounds: {}", src.model.bounds_hi.transpose());

            auto fit_filter =
                QI::ModelFitFilter<FitType>::New(&src, verbose, covar, resids, subregion.Get());
            fit_filter->ReadInputs(
                {spgr_path.Get(), ssfp_path.Get()}, {f0.Get(), B1.Get()}, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + model_name);
            QI::Log(verbose, "Finished.");
        }
    };
    switch (modelarg.Get()) {
    case 2: {
        QI::TwoPoolModel two_pool{spgr, ssfp, scale.Get()};
        process(two_pool, "2C_");
    } break;
    case 3: {
        QI::ThreePoolModel three_pool{spgr, ssfp, scale.Get()};
        process(three_pool, "3C_");
    } break;
    default:
        QI::Fail("Unknown model specifier: {}", modelarg.Get());
    }
    return EXIT_SUCCESS;
}
