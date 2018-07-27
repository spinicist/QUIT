/*
 *  qi_ase_oef.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>

#include "Model.h"
#include "FitFunction.h"
#include "ModelFitFilter.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "MultiEchoSequence.h"

using namespace std::literals;

using ASEModel = QI::Model<4, 3, QI::MultiEchoFlexSequence>;
template<> std::array<const std::string, 4> ASEModel::varying_names{{"R2prime"s, "DBV"s, "OEF"s, "dHb"s}};
template<> std::array<const std::string, 3> ASEModel::fixed_names{{"grad_x"s, "grad_y"s, "grad_z"s}};
template<> const QI_ARRAYN(double, 3) ASEModel::fixed_defaults{0.0, 0.0, 0.0};

using ASEFit = QI::FitFunction<ASEModel>;

struct ASELLS : ASEFit {
    const double B0;
    const QI::VolumeF::SpacingType voxsize;
    Eigen::ArrayXd TE_above_Tc;
    Eigen::ArrayXi TE_indices, TE0_indices;
    // Constants for calculations
    const double kappa = 0.03; // Conversion factor
    const double gamma = 42.577e6; // Gyromagnetic Ratio
    const double delta_X0 = 0.264e-6; // Difference in susceptibility of oxy and fully de-oxy blood
    const double Hb = 0.34 / kappa; // Hct = 0.34;

    ASELLS(QI::MultiEchoFlexSequence *seq, const double B0, const QI::VolumeF::SpacingType voxsize) :
        ASEFit(seq), B0(B0), voxsize(voxsize)
    {
        // Nic Blockley uses Tc = 15 ms for 3T, scale for other field-strengths
        const int Tc = 0.015 / (B0 / 3);
        const int above_Tc_count = (seq->TE.abs() > Tc).count();
        const int zero_count = (seq->TE.abs() == 0.0).count();
        if (zero_count == 0) {
            QI_FAIL("Did not find a zero echo-time in input");
        }
        if (above_Tc_count == 0) {
            QI_FAIL("No echo-times above critical value");
        }
        TE_above_Tc = Eigen::ArrayXd(above_Tc_count);
        TE_indices = Eigen::ArrayXi(above_Tc_count);
        TE0_indices = Eigen::ArrayXi(zero_count);
        int zero_index = 0, index = 0;
        for (Eigen::Index i = 0; i < seq->size(); i++) {
            if (seq->TE(i) == 0.0) {
                TE0_indices(zero_index) = i;
                zero_index++;
            }
            if (std::abs(seq->TE(i)) > Tc) {
                TE_indices(index) = i;
                TE_above_Tc(index) = std::abs(seq->TE(i));
                index++;
            }
        }
    }


    double sinc(const double x) const {
        static double const taylor_0_bound = std::numeric_limits<double>::epsilon();
        static double const taylor_2_bound = sqrt(taylor_0_bound);
        static double const taylor_n_bound = sqrt(taylor_2_bound);

        if (std::abs(x) >= taylor_n_bound) {
            return(sin(x)/x);
        } else {
            // approximation by taylor series in x at 0 up to order 0
            double result = 1;

            if (abs(x) >= taylor_0_bound) {
                double x2 = x*x;
                // approximation by taylor series in x at 0 up to order 2
                result -= x2/6;

                if (abs(x) >= taylor_2_bound) {
                    // approximation by taylor series in x at 0 up to order 4
                    result += (x2*x2)/120;
                }
            }
            return result;
        }
    }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, ASEModel::NV) &outputs,
                          ResidualType &/*Unused*/, std::vector<Eigen::ArrayXd> &/*Unused*/, FlagType &/*Unused*/) const override
    {
        const Eigen::ArrayXd &all_data = inputs[0];
        Eigen::ArrayXd data = Eigen::ArrayXd(TE_indices.rows());
        for (Eigen::Index i = 0; i < TE_indices.rows(); i++) {
            double F = 1.0;
            for (auto d = 0; d < 3; d++) {
                const double grad = fixed[d];
                const double x = grad * 2 * M_PI * voxsize[d] * sequence->TE[TE_indices[i]] / 2;
                F *= std::abs(sinc(x));
            }
            data[i] = all_data[TE_indices[i]] / F;
        }

        Eigen::MatrixXd X(TE_above_Tc.rows(), 2);
        X.col(0) = TE_above_Tc;
        X.col(1).setOnes();
        Eigen::VectorXd Y = data.log();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        double sumTE0 = 0;
        for (int i = 0; i < TE0_indices.rows(); i++) {
            sumTE0 += data(TE0_indices[i]);
        }
        const double logTE0 = log(sumTE0 / TE0_indices.rows());
        const double R2prime = -b[0];
        const double logS0_linear = b[1];
        const double DBV = logS0_linear - logTE0;
        const double dHb = 3*R2prime / (DBV * 4 * gamma * M_PI * delta_X0 * kappa * B0);
        const double OEF = dHb / Hb;

        outputs[0] = R2prime;
        outputs[1] = DBV*100;
        outputs[2] = OEF*100;
        outputs[3] = dHb;
        return std::make_tuple(true, "");
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates the OEF from ASE data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ASE_FILE", "Input ASE file");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<double> B0(parser, "B0", "Field-strength (Tesla), default 3", {'B', "B0"}, 3.0);
    args::ValueFlag<std::string> gradx(parser, "GRADX", "Gradient of field-map in x-direction for MFG correction", {'x', "gradx"});
    args::ValueFlag<std::string> grady(parser, "GRADY", "Gradient of field-map in y-direction for MFG correction", {'y', "grady"});
    args::ValueFlag<std::string> gradz(parser, "GRADZ", "Gradient of field-map in z-direction for MFG correction", {'z', "gradz"});
    args::ValueFlag<double> slice_arg(parser, "SLICE THICKNESS", "Slice-thickness for MFG calculation (useful if there was a slice gap)", {'s', "slice"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI_LOG(verbose, "Reading ASE data from: " << QI::CheckPos(input_path));
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path));
    rapidjson::Document json = QI::ReadJSON(std::cin);
    QI::MultiEchoFlexSequence sequence(json["MultiEchoFlex"]);
    QI::VolumeF::SpacingType vox_size = input->GetSpacing();
    if (slice_arg) {
        vox_size[2]  = slice_arg.Get();
    }
    ASELLS fit(&sequence, B0.Get(), vox_size);
    auto fit_filter = itk::ModelFitFilter<ASEFit>::New();
    fit_filter->SetVerbose(verbose);
    fit_filter->SetFitFunction(&fit);
    fit_filter->SetOutputAllResiduals(false);
    QI_LOG(verbose, "Using " << threads.Get() << " threads" );
    fit_filter->SetInput(0, input);
    if (mask) fit_filter->SetMask(QI::ReadImage(mask.Get()));
    if (gradx) fit_filter->SetFixed(0, QI::ReadImage(gradx.Get()));
    if (grady) fit_filter->SetFixed(1, QI::ReadImage(grady.Get()));
    if (gradz) fit_filter->SetFixed(2, QI::ReadImage(gradz.Get()));
    if (subregion) {
        fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
    }
    QI_LOG(verbose, "Processing");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        fit_filter->AddObserver(itk::ProgressEvent(), monitor);
    }
    fit_filter->Update();
    QI_LOG(verbose, "Elapsed time was " << fit_filter->GetTotalTime() << "s");
    
    for (size_t i = 0; i < ASEModel::NV; i++) {
        const std::string fname = outPrefix + "_" + ASEModel::varying_names[i] + QI::OutExt();
        QI_LOG(verbose, "Writing file: " << fname);
        QI::WriteImage(fit_filter->GetOutput(i), fname);
    }
    return EXIT_SUCCESS;
}
