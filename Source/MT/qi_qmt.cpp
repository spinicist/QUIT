/*
 *  qi_qmt.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood, Samuel Hurley, Erika Raven
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
#include "FitScaledAuto.h"
#include "ImageIO.h"
#include "Lineshape.h"
#include "MTSequences.h"
#include "Macro.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct RamaniModel : QI::Model<double, double, 5, 3, 1, 3> {
    QI::ZSpecSequence const &                  sequence;
    const QI::Lineshapes                       lineshape;
    const std::shared_ptr<QI::InterpLineshape> interp = nullptr;

    std::array<const std::string, NV> const varying_names{
        {"M0_f"s, "F_over_R1_f"s, "T2_b"s, "T1_f_over_T2_f"s, "k"s}};
    std::array<const std::string, 3> const  derived_names{{"T1_f"s, "T2_f"s, "f_b"s}};
    std::array<const std::string, NF> const fixed_names{{"f0"s, "B1"s, "T1_app"}};

    FixedArray const   fixed_defaults{0.0, 1.0, 1.0};
    VaryingArray const bounds_lo{0.5, 1.e-12, 1.e-7, 1, 1};
    VaryingArray const bounds_hi{1.5, 1, 100.e-6, 100., 1e2};
    VaryingArray const start{1., 0.1, 10.e-6, 25., 10.};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v, const FixedArray &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        // Don't use Ramani's notation
        const auto &R1_b           = 2.5;  // Fix
        const auto &M0_f           = v[0]; // We normalise out the gain in the fit-function
        const auto &F_over_R1_f    = v[1];
        const auto &T2b            = v[2];
        const auto &T1_f_over_T2_f = v[3];
        const auto &k              = v[4];
        const auto &f0             = f[0];
        const auto &B1             = f[1];
        const auto  kM0_f          = k * M0_f;
        QI_ARRAY(typename Derived::Scalar) lsv;
        switch (lineshape) {
        case QI::Lineshapes::Gaussian:
            lsv = QI::Gaussian((sequence.sat_f0 + f0), T2b);
            break;
        case QI::Lineshapes::Lorentzian:
            lsv = QI::Lorentzian((sequence.sat_f0 + f0), T2b);
            break;
        case QI::Lineshapes::SuperLorentzian:
            lsv = QI::SuperLorentzian((sequence.sat_f0 + f0), T2b);
            break;
        case QI::Lineshapes::Interpolated:
            lsv = (*interp)((sequence.sat_f0 + f0), T2b);
            break;
        }

        const auto w_cwpe = (B1 * sequence.sat_angle / sequence.pulse.p1) *
                            sqrt(sequence.pulse.p2 / (sequence.Trf * sequence.TR));
        const auto R_rfb = M_PI * (w_cwpe * w_cwpe) * lsv;

        const auto S = M0_f * (R1_b * kM0_f * F_over_R1_f + R_rfb + R1_b + kM0_f) /
                       ((kM0_f * F_over_R1_f) * (R1_b + R_rfb) +
                        (1.0 + pow(w_cwpe / (2 * M_PI * sequence.sat_f0), 2.0) * T1_f_over_T2_f) *
                            (R_rfb + R1_b + kM0_f));
        QI_DBVEC(v)
        QI_DBVEC(w_cwpe)
        QI_DBVEC(R_rfb)
        QI_DBVEC(S)

        return S;
    }

    void derived(const VaryingArray &p, const FixedArray &f, DerivedArray &d) const {
        // Convert from the fitted parameters to useful ones
        const auto &R1_b = 2.5; // Fix
        // const auto &M0_f           = p[0];
        const auto &F_over_R1_f = p[1];
        // const auto &T2b            = p[2];
        const auto &T1_f_over_T2_f = p[3];
        const auto &k              = p[4];
        const auto &T1_obs         = f[2];

        const auto R_obs = 1 / T1_obs;
        // ((R1_b < R_obs) && (R1_b - R_obs + k) > 0) ?
        //                     R_obs :
        const auto Ra  = R_obs / (1.0 + ((k * F_over_R1_f * (R1_b - R_obs)) / (R1_b - R_obs + k)));
        const auto F   = F_over_R1_f * Ra;
        const auto f_b = F / (1.0 - F);

        //{"T1_f"s, "T2_f"s, "f_b"s}
        d[0] = QI::Clamp(1.0 / Ra, 0., 5.0);
        d[1] = QI::Clamp(F_over_R1_f / (F * T1_f_over_T2_f), 0., 3.);
        d[2] = QI::Clamp(100. * f_b, 0., 100.);
        QI_DBVEC(p)
        QI_DB(R_obs)
        QI_DB(Ra)
        QI_DBVEC(d)
    }
};

using RamaniFitFunction = QI::ScaledAutoDiffFit<RamaniModel>;

//******************************************************************************
// Main
//******************************************************************************
int qmt_main(args::Subparser &parser) {
    args::Positional<std::string> mtsat_path(parser, "MTSAT FILE", "Path to MT-Sat data");
    QI_COMMON_ARGS;
    args::ValueFlag<std::string> T1(parser, "T1", "T1 map (seconds) file ** REQUIRED **", {"T1"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (Hz) file", {'f', "f0"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> lineshape_arg(
        parser,
        "LINESHAPE",
        "Either Gaussian, Lorentzian, Superlorentzian, or a .json file generated by qi_lineshape",
        {'l', "lineshape"},
        "Gaussian");
    parser.Parse();
    QI::CheckPos(mtsat_path);
    QI::Log(verbose, "Reading sequence information");
    json           input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto           mtsat_sequence = input.at("MTSat").get<QI::ZSpecSequence>();
    QI::Lineshapes lineshape;
    std::shared_ptr<QI::InterpLineshape> interp = nullptr;
    if (lineshape_arg.Get() == "Gaussian") {
        QI::Log(verbose, "Using a Gaussian lineshape");
        lineshape = QI::Lineshapes::Gaussian;
    } else if (lineshape_arg.Get() == "Lorentzian") {
        QI::Log(verbose, "Using a Lorentzian lineshape");
        lineshape = QI::Lineshapes::Lorentzian;
    } else if (lineshape_arg.Get() == "Superlorentzian") {
        QI::Log(verbose, "Using a Super-Lorentzian lineshape");
        lineshape = QI::Lineshapes::SuperLorentzian;
    } else {
        QI::Log(verbose, "Reading lineshape file: {}", lineshape_arg.Get());
        json ls_file = QI::ReadJSON(lineshape_arg.Get());
        interp       = std::make_shared<QI::InterpLineshape>(
            ls_file.at("lineshape").get<QI::InterpLineshape>());
        lineshape = QI::Lineshapes::Interpolated;
    }

    RamaniModel model{{}, mtsat_sequence, lineshape, interp};
    if (simulate) {
        QI::SimulateModel<RamaniModel, false>(input,
                                              model,
                                              {f0.Get(), B1.Get(), QI::CheckValue(T1)},
                                              {mtsat_path.Get()},
                                              mask.Get(),
                                              verbose,
                                              simulate.Get(),
                                              subregion.Get());
    } else {
        RamaniFitFunction fit{model};

        auto fit_filter = QI::ModelFitFilter<RamaniFitFunction>::New(
            &fit, verbose, covar, resids, subregion.Get());
        fit_filter->ReadInputs(
            {mtsat_path.Get()}, {f0.Get(), B1.Get(), QI::CheckValue(T1)}, mask.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "QMT_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
