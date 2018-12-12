/*
 *  qi_asl.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <iostream>

#include "Args.h"
#include "CASLSequence.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"
#include "itkUnaryFunctorImageFilter.h"

using namespace std::literals;

using CASLModel = QI::Model<1, 2, QI::CASLSequence>;
template <> std::array<const std::string, 1> CASLModel::varying_names{{"CBF"s}};
template <>
std::array<const std::string, 2> CASLModel::fixed_names{{"T1_tisse"s, "PD"s}};
template <> const QI_ARRAYN(double, 2) CASLModel::fixed_defaults{0.0, 0.0};

struct CASLFit {
    static const bool Blocked = true;
    static const bool Indexed = true;
    using InputType           = double;
    using OutputType          = double;
    using ResidualType        = double;
    using FlagType            = int; // Iterations
    using ModelType           = CASLModel;
    using SequenceType        = QI::CASLSequence;
    ModelType &  model;
    const double T1_blood, alpha, lambda;

    CASLFit(ModelType &c, const double &t, const double &a, const double &l)
        : model(c), T1_blood(t), alpha(a), lambda(l) {}
    int n_inputs() const { return 1; }
    int input_size(const int /* Unused */) const { return 2; }
    int n_fixed() const { return ModelType::NF; }
    int n_outputs() const { return ModelType::NV; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, CASLModel::NV) & outputs,
                          ResidualType & /*Unused*/,
                          std::vector<Eigen::ArrayXd> & /*Unused*/,
                          FlagType & /*Unused*/, const int /*Unused*/,
                          const itk::Index<3> &voxel) const {
        const double label   = inputs[0][0];
        const double control = inputs[0][1];

        const double diff      = control - label;
        const double T1_tissue = fixed[0];
        const double PD_correction =
            (T1_tissue == 0.0) ? 1.0
                               : (1 - exp(-model.sequence.TR / T1_tissue));
        const double PD =
            ((fixed[1] == 0.0) ? control : fixed[1]) / PD_correction;
        const double PLD = (model.sequence.post_label_delay.rows() > 1)
                               ? model.sequence.post_label_delay[voxel[2]]
                               : model.sequence.post_label_delay[0];
        const double CBF = (6000 * lambda * diff * exp(PLD / T1_blood)) /
                           (2. * alpha * T1_blood * PD *
                            (1. - exp(-model.sequence.label_time / T1_blood)));
        outputs[0] = CBF;
        return std::make_tuple(true, "");
    }
};

struct VectorMean {
    VectorMean(){};
    ~VectorMean(){};
    bool operator!=(const VectorMean &) const { return false; }
    bool operator==(const VectorMean &other) const { return !(*this != other); }
    inline float operator()(const itk::VariableLengthVector<float> &vec) const {
        float sum = 0;
        for (size_t i = 0; i < vec.Size(); i++) {
            sum += vec[i];
        }
        return sum / vec.Size();
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates CBF from ASL data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ASL_FILE",
                                             "Input ASL file");
    args::HelpFlag                help(parser, "HELP", "Show this help message",
                        {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information",
                       {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION",
        "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::Flag average(parser, "AVERAGE", "Average the time-series",
                       {'a', "average"});
    args::Flag slice_time(parser, "SLICE TIME CORRECTION",
                          "Apply slice-time correction (number of post-label "
                          "delays must match number of slices)",
                          {'s', "slicetime"});
    args::ValueFlag<double> T1_blood(
        parser, "BLOOD T1",
        "Value of blood T1 to use (seconds), default 1.65 for 3T",
        {'b', "blood"}, 1.65);
    args::ValueFlag<std::string> T1_tissue_path(
        parser, "TISSUE T1", "Path to tissue T1 map (units are seconds)",
        {'t', "tissue"});
    args::ValueFlag<std::string> PD_path(parser, "PROTON DENSITY",
                                         "Path to PD image", {'p', "pd"});
    args::ValueFlag<double>      alpha(parser, "ALPHA",
                                  "Labelling efficiency, default 0.9",
                                  {'a', "alpha"}, 0.9);
    args::ValueFlag<double>      lambda(
        parser, "LAMBDA", "Blood-brain partition co-efficent, default 0.9 mL/g",
        {'l', "lambda"}, 0.9);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path), verbose);
    QI_LOG(verbose, "Reading sequence parameters");
    rapidjson::Document json = QI::ReadJSON(std::cin);
    QI::CASLSequence    sequence(json["CASL"]);
    CASLModel           model{sequence};
    CASLFit casl_fit{model, T1_blood.Get(), alpha.Get(), lambda.Get()};
    auto    fit_filter = itk::ModelFitFilter<CASLFit>::New(&casl_fit);
    fit_filter->SetVerbose(verbose);
    fit_filter->SetInput(0, input);
    if (PD_path)
        fit_filter->SetFixed(0, QI::ReadImage(PD_path.Get(), verbose));
    if (T1_tissue_path)
        fit_filter->SetFixed(1, QI::ReadImage(T1_tissue_path.Get(), verbose));
    const auto n_slices = input->GetLargestPossibleRegion().GetSize()[2];
    if (slice_time &&
        (n_slices != static_cast<size_t>(sequence.post_label_delay.rows()))) {
        QI_FAIL("Number of post-label delays "
                << sequence.post_label_delay.rows()
                << " does not match number of slices " << n_slices);
    }
    fit_filter->SetBlocks(input->GetNumberOfComponentsPerPixel() / 2);
    if (mask)
        fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
    if (subregion)
        fit_filter->SetSubregion(QI::RegionArg(subregion.Get()));
    QI_LOG(verbose, "Processing");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        fit_filter->AddObserver(itk::ProgressEvent(), monitor);
    }
    fit_filter->Update();
    QI_LOG(verbose, "Elapsed time was " << fit_filter->GetTotalTime() << "s\n"
                                        << "Writing results files.");
    std::string outPrefix = outarg.Get() + "CASL";
    if (average) {
        auto mean_filter =
            itk::UnaryFunctorImageFilter<QI::VectorVolumeF, QI::VolumeF,
                                         VectorMean>::New();
        mean_filter->SetInput(fit_filter->GetOutput(0));
        mean_filter->Update();
        QI::WriteImage(mean_filter->GetOutput(),
                       outPrefix + "_CBF" + QI::OutExt());
    } else {
        QI::WriteVectorImage(fit_filter->GetOutput(0),
                             outPrefix + "_CBF" + QI::OutExt());
    }
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
