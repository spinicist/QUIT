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

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "itkTimeProbe.h"

#include "Types.h"
#include "Util.h"
#include "Args.h"
#include "IO.h"
#include "Model.h"
#include "SequenceGroup.h"
#include "RegionContraction.h"
#include "ApplyAlgorithmFilter.h"
#include "ReorderVectorFilter.h"

struct MCDSRCFunctor {
    const QI::SequenceGroup &m_sequence;
    const Eigen::ArrayXd m_data, m_weights;
    const std::shared_ptr<QI::Model> m_model;

    MCDSRCFunctor(std::shared_ptr<QI::Model> m, QI::SequenceGroup &s,
                  const Eigen::ArrayXd &d, const Eigen::ArrayXd &w) :
        m_sequence(s), m_data(d), m_model(m), m_weights(w)
    {
        assert(static_cast<size_t>(m_data.rows()) == m_sequence.size());
    }

    int inputs() const { return m_model->nParameters(); }
    int values() const { return m_sequence.size(); }

    const bool constraint(const Eigen::VectorXd &params) const {
        return m_model->ValidParameters(params);
    }

    Eigen::ArrayXd residuals(const Eigen::Ref<Eigen::VectorXd> &params) const {
        const Eigen::ArrayXd s = (m_sequence.signal(m_model, params)).abs();
        return m_data - s;
    }
    double operator()(const Eigen::Ref<Eigen::VectorXd> &params) const {
        return (residuals(params) * m_weights).square().sum();
    }
};

struct SRCAlgo : public QI::ApplyF::Algorithm {
    Eigen::ArrayXXd m_bounds;
    std::shared_ptr<QI::Model> m_model = nullptr;
    QI::SequenceGroup &m_sequence;
    QI::FieldStrength m_tesla = QI::FieldStrength::Three;
    int m_iterations = 0;
    size_t m_samples = 5000, m_retain = 50;
    bool m_gauss = true;

    SRCAlgo(std::shared_ptr<QI::Model>&m, Eigen::ArrayXXd &b,
            QI::SequenceGroup &s, int mi) :
        m_model(m), m_bounds(b), m_sequence(s), m_iterations(mi)
    {}

    size_t numInputs() const override  { return m_sequence.count(); }
    size_t numOutputs() const override { return m_model->nParameters(); }
    size_t dataSize() const override   { return m_sequence.size(); }

    void setModel(std::shared_ptr<QI::Model> &m) { m_model = m; }
    void setSequence(QI::SequenceGroup &s) { m_sequence = s; }
    void setBounds(Eigen::ArrayXXd &b) { m_bounds = b; }
    void setIterations(const int i) { m_iterations = i; }
    float zero() const override { return 0.f; }

    void setGauss(bool g) { m_gauss = g; }
    size_t numConsts() const override  { return 2; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(2);
        def[0] = NAN; def[1] = 1.0f; // f0, B1
        return def;
    }

    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIters &its) const override
    {
        Eigen::ArrayXd data(dataSize());
        int dataIndex = 0;
        for (int i = 0; i < inputs.size(); i++) {
            Eigen::Map<const Eigen::ArrayXf> this_data(inputs[i].GetDataPointer(), inputs[i].Size());
            if (m_model->scaleToMean()) {
                data.segment(dataIndex, this_data.rows()) = this_data.cast<double>() / this_data.abs().mean();
            } else {
                data.segment(dataIndex, this_data.rows()) = this_data.cast<double>();
            }
            dataIndex += this_data.rows();
        }
        Eigen::ArrayXd thresh(m_model->nParameters()); thresh.setConstant(0.05);
        const double f0 = consts[0];
        const double B1 = consts[1];
        Eigen::ArrayXXd localBounds = m_bounds;
        Eigen::ArrayXd weights = Eigen::ArrayXd::Ones(m_sequence.size());
        if (isfinite(f0)) { // We have an f0 map, add it to the fitting bounds
            localBounds.row(m_model->ParameterIndex("f0")) += f0;
            weights = m_sequence.weights(f0);
        }
        localBounds.row(m_model->ParameterIndex("B1")).setConstant(B1);
        MCDSRCFunctor func(m_model, m_sequence, data, weights);
        QI::RegionContraction<MCDSRCFunctor> rc(func, localBounds, thresh, m_samples, m_retain, m_iterations, 0.02, m_gauss, false);
        Eigen::ArrayXd pars(m_model->nParameters());
        rc.optimise(pars);
        for (int i = 0; i < m_model->nParameters(); i++) {
            outputs[i] = pars[i];
        }
        Eigen::ArrayXf r = func.residuals(pars).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = rc.contractions();
        return true;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates MWF & other parameter maps from mcDESPOT data\n"
                                "All times (e.g. T1, TR) are in SECONDS. All angles are in degrees.\n"
                                "http://github.com/spinicist/QUIT");
    args::PositionalList<std::string> input_paths(parser, "INPUT FILES", "Input image files");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (Hertz)", {'f', "f0"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<std::string> modelarg(parser, "MODEL", "Select model to fit - 1/2/2nex/3/3_f0/3nex, default 3", {'M', "model"}, "3");
    args::Flag scale(parser, "SCALE", "Normalize signals to mean (a good idea)", {'S', "scale"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Select (S)tochastic or (G)aussian Region Contraction", {'a', "algo"}, 'G');
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations, default 4", {'i',"its"}, 4);
    args::ValueFlag<char> field(parser, "FIELD STRENGTH", "Specify field-strength for fitting regions - 3/7/u for user input", {'t', "tesla"}, '3');
    QI::ParseArgs(parser, argc, argv);

    std::shared_ptr<QI::Model> model = nullptr;
    if (modelarg.Get() == "1")         { model = std::make_shared<QI::SCD>(); }
    else if (modelarg.Get() == "2")    { model = std::make_shared<QI::MCD2>(); }
    else if (modelarg.Get() == "2nex") { model = std::make_shared<QI::MCD2_NoEx>(); }
    else if (modelarg.Get() == "3")    { model = std::make_shared<QI::MCD3>(); }
    else if (modelarg.Get() == "3_f0") { model = std::make_shared<QI::MCD3_f0>(); }
    else if (modelarg.Get() == "3nex") { model = std::make_shared<QI::MCD3_NoEx>(); }
    else {
        std::cerr << "Invalid model " << modelarg.Get() << " specified." << std::endl;
        return EXIT_FAILURE;
    }
    if (verbose) std::cout << "Using " << model->Name() << " model." << std::endl;
    if (verbose && scale) std::cout << "Mean-scaling selected" << std::endl;
    model->setScaleToMean(scale);

    Eigen::ArrayXXd bounds;
    Eigen::ArrayXd start;
    switch (field.Get()) {
        case '3':
            bounds = model->Bounds(QI::FieldStrength::Three);
            start = model->Default(QI::FieldStrength::Three);
            break;
        case '7':
            bounds = model->Bounds(QI::FieldStrength::Seven);
            start = model->Bounds(QI::FieldStrength::Seven);
            break;
        case 'u': {
            Eigen::ArrayXd temp;
            if (verbose) std::cout << "Enter lower bounds" << std::endl;
            QI::ReadArray(std::cin, temp);
            bounds.col(0) = temp;
            if (verbose) std::cout << "Enter upper bounds" << std::endl;
            QI::ReadArray(std::cin, temp);
            bounds.col(1) = temp;
        } break;
    default:
        std::cerr << "Unknown boundaries type " << field.Get() << std::endl;
        return EXIT_FAILURE;
        break;
    }

    std::vector<QI::VectorVolumeF::Pointer> images(input_paths.Get().size());
    for (auto &input_path : input_paths.Get()) {
        if (verbose) std::cout << "Reading file: " << input_path << std::endl;
        auto image = QI::ReadVectorImage<float>(input_path);
        image->DisconnectPipeline(); // This step is really important.
        images.push_back(image);
    }
    auto sequences = QI::ReadSequence<QI::SequenceGroup>(std::cin, "SequenceGroup", verbose);
    if (sequences.count() != images.size()) {
        QI_FAIL("Sequence group size " << sequences.count() << " does not match images size " << images.size());
    }
    auto apply = QI::ApplyF::New();
    switch (algorithm.Get()) {
        case 'S': {
            if (verbose) std::cout << "Using SRC algorithm" << std::endl;
            std::shared_ptr<SRCAlgo> algo = std::make_shared<SRCAlgo>(model, bounds, sequences, its.Get());
            algo->setGauss(false);
            apply->SetAlgorithm(algo);
        } break;
        case 'G': {
            if (verbose) std::cout << "Using GRC algorithm" << std::endl;
            std::shared_ptr<SRCAlgo> algo = std::make_shared<SRCAlgo>(model, bounds, sequences, its.Get());
            algo->setGauss(true);
            apply->SetAlgorithm(algo);
        } break;
        default:
            std::cerr << "Unknown algorithm type " << algorithm.Get() << std::endl;
            return EXIT_FAILURE;
    }
    apply->SetOutputAllResiduals(resids);
    apply->SetVerbose(verbose);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()); // mcdespot with a mask & threads is a very unbalanced algorithm
    for (int i = 0; i < images.size(); i++) {
        apply->SetInput(i, images[i]);
    }
    if (f0) apply->SetConst(0, QI::ReadImage(f0.Get()));
    if (B1) apply->SetConst(1, QI::ReadImage(B1.Get()));
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (subregion) apply->SetSubregion(QI::RegionArg(args::get(subregion)));

    // Need this here so the bounds.txt file will have the correct prefix
    std::string outPrefix = outarg.Get() + model->Name() + "_";
    if (verbose) {
        std::cout << "Bounds:\n" <<  bounds.transpose() << std::endl;
        std::ofstream boundsFile(outPrefix + "bounds.txt");
        boundsFile << "Names: ";
        for (size_t p = 0; p < model->nParameters(); p++) {
            boundsFile << model->ParameterNames()[p] << "\t";
        }
        boundsFile << std::endl << "Bounds:\n" << bounds.transpose() << std::endl;
        boundsFile.close();
    }

    if (verbose) {
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;
    }
    for (int i = 0; i < model->nParameters(); i++) {
        QI::WriteImage(apply->GetOutput(i), outPrefix + model->ParameterNames()[i] + QI::OutExt());
    }
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual" + QI::OutExt());
    if (resids) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals" + QI::OutExt());
    }
    QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "iterations" + QI::OutExt());
    return EXIT_SUCCESS;
}

