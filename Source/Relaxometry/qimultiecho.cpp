/*
 *  multiecho.cpp
 *
 *  Created by Tobias Wood on 27/01/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "itkTimeProbe.h"
#include "itkImageFileReader.h"
#include "itkTileImageFilter.h"

#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "MultiEchoSequence.h"
#include "ApplyTypes.h"
#include "ImageToVectorFilter.h"

typedef itk::ImageToVectorFilter<QI::SeriesF> SeriesToVectorF;

/*
 * Base class for the 3 different algorithms
 */
class RelaxAlgo : public QI::ApplyF::Algorithm {
private:
    const std::shared_ptr<QI::Model::OnePool> m_model = std::make_shared<QI::Model::OnePool>();
protected:
    QI::MultiEchoSequence m_sequence;
    double m_clampLo = -std::numeric_limits<double>::infinity();
    double m_clampHi = std::numeric_limits<double>::infinity();
    double m_thresh = -std::numeric_limits<double>::infinity();
    size_t m_iterations = 15;

    void clamp_and_threshold(const Eigen::ArrayXd &data, std::vector<TOutput> &outputs, TConst &residual, TInput &resids,
                            const double PD, const double T2) const {
        if (PD > m_thresh) {
            outputs[0] = PD;
            outputs[1] = QI::Clamp(T2, m_clampLo, m_clampHi);
            Eigen::ArrayXd theory = QI::One_MultiEcho(m_sequence.TE, m_sequence.TR, PD, 0., T2).array().abs(); // T1 isn't modelled, set to 0 for instant recovery
            Eigen::ArrayXf r = (data.array() - theory).cast<float>();
            residual = sqrt(r.square().sum() / r.rows());
            resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        } else {
            outputs[0] = 0.;
            outputs[1] = 0.;
            resids = itk::VariableLengthVector<float>(data.rows());
            resids.Fill(0);
        }
    }
public:
    void setIterations(size_t n) { m_iterations = n; }
    void setSequence(const QI::MultiEchoSequence &s) { m_sequence = s; }
    void setClamp(double lo, double hi) { m_clampLo = lo; m_clampHi = hi; }
    void setThresh(double t) { m_thresh = t; }
    size_t numInputs() const override { return m_sequence.count(); }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 2; }
    size_t dataSize() const override { return m_sequence.size(); }
    float zero() const override { return 0.f; }

    std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 1.0); // B1
        return def;
    }
};

class LogLinAlgo: public RelaxAlgo {
public:
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> & /* Unused */,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        // Set up echo times array
        Eigen::MatrixXd X(m_sequence.size(), 2);
        X.col(0) = m_sequence.TE;
        X.col(1).setOnes();
        Eigen::VectorXd Y = data.array().log();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        double PD = exp(b[1]);
        double T2 = -1 / b[0];
        clamp_and_threshold(data, outputs, residual, resids, PD, T2);
        its = 1;
        return std::make_tuple(true, "");
    }
};

class ARLOAlgo : public RelaxAlgo {
public:
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> & /* Unused */,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        const double dTE_3 = (m_sequence.ESP / 3);
        double si2sum = 0, di2sum = 0, sidisum = 0;
        for (Eigen::Index i = 0; i < m_sequence.size() - 2; i++) {
            const double si = dTE_3 * (data(i) + 4*data(i+1) + data(i+2));
            const double di = data(i) - data(i+2);
            si2sum += si*si;
            di2sum += di*di;
            sidisum += si*di;
        }
        double T2 = (si2sum + dTE_3*sidisum) / (dTE_3*di2sum + sidisum);
        double PD = (data.array() / exp(-m_sequence.TE / T2)).mean();
        clamp_and_threshold(data, outputs, residual, resids, PD, T2);
        its = 1;
        return std::make_tuple(true, "");
    }
};

class RelaxFunctor : public Eigen::DenseFunctor<double> {
    protected:
        const QI::MultiEchoSequence m_sequence;
        const Eigen::ArrayXd m_data;
        const std::shared_ptr<QI::Model::OnePool> m_model = std::make_shared<QI::Model::OnePool>();

    public:
        RelaxFunctor(const QI::MultiEchoSequence &cs, const Eigen::ArrayXd &data) :
            DenseFunctor<double>(2, cs.size()),
            m_sequence(cs), m_data(data)
        {
            assert(static_cast<size_t>(m_data.rows()) == values());
        }

        int operator()(const Eigen::Ref<Eigen::VectorXd> &params, Eigen::Ref<Eigen::ArrayXd> diffs) const {
            eigen_assert(diffs.size() == values());
            Eigen::VectorXd fullp(5);
            fullp << params(0), 0, params(1), 0, 1.0; // Fix B1 to 1.0 for now
            Eigen::ArrayXcd s = m_sequence.signal(m_model, fullp);
            diffs = s.abs() - m_data;
            return 0;
        }
};

class NonLinAlgo : public RelaxAlgo {
public:
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> & /* Unused */,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        const Eigen::ArrayXd data = indata.cast<double>();
        RelaxFunctor f(m_sequence, data);
        Eigen::NumericalDiff<RelaxFunctor> nDiff(f);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<RelaxFunctor>> lm(nDiff);
        lm.setMaxfev(m_iterations * (m_sequence.size() + 1));
        // Just PD & T2 for now
        // Basic guess of T2=50ms
        Eigen::VectorXd p(2); p << data(0), 0.05;
        lm.minimize(p);
        clamp_and_threshold(data, outputs, residual, resids, p[0], p[1]);
        its = lm.iterations();
        return std::make_tuple(true, "");
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates T2/T2* maps from multi-echo data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT FILE", "Input multi-echo data");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/a/n)", {'a',"algo"}, 'l');
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations for WLLS/NLLS (default 15)", {'i',"its"}, 15);
    args::ValueFlag<float> clampT2(parser, "CLAMP T2", "Clamp T2 between 0 and value", {'p',"clampPD"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<float> threshPD(parser, "THRESHOLD PD", "Only output maps when PD exceeds threshold value", {'t', "tresh"});
    QI::ParseArgs(parser, argc, argv, verbose);

    QI_LOG(verbose, "Opening input file: " << QI::CheckPos(input_path));
    auto inputFile = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path));

    std::shared_ptr<RelaxAlgo> algo = ITK_NULLPTR;
    switch (algorithm.Get()) {
        case 'l': algo = std::make_shared<LogLinAlgo>(); QI_LOG(verbose, "LogLin algorithm selected." ); break;
        case 'a': algo = std::make_shared<ARLOAlgo>(); QI_LOG(verbose, "ARLO algorithm selected." ); break;
        case 'n': algo = std::make_shared<NonLinAlgo>(); QI_LOG(verbose, "Non-linear algorithm (Levenberg Marquardt) selected." ); break;
        default:
            QI_FAIL("Unknown algorithm type " << algorithm.Get());
    }
    algo->setThresh(threshPD.Get());
    algo->setClamp(0, clampT2.Get());
    algo->setIterations(its.Get());

    // Gather input data
    QI_LOG(verbose, "Reading sequence parameters");
    rapidjson::Document input = QI::ReadJSON(std::cin);
    QI::MultiEchoSequence multiecho(input["MultiEcho"]);
    size_t nVols = inputFile->GetLargestPossibleRegion().GetSize()[3] / multiecho.size();
    algo->setSequence(multiecho);
    auto apply = QI::ApplyF::New();
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()); // Unbalanced algorithm
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));

    auto PDoutput = itk::TileImageFilter<QI::VolumeF, QI::SeriesF>::New();
    auto T2output = itk::TileImageFilter<QI::VolumeF, QI::SeriesF>::New();
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1; layout[3] = nVols;
    PDoutput->SetLayout(layout);
    T2output->SetLayout(layout);
    QI_LOG(verbose, "Processing" );
    auto inputVector = SeriesToVectorF::New();
    inputVector->SetInput(inputFile);
    inputVector->SetBlockSize(multiecho.size());
    std::vector<QI::VolumeF::Pointer> PDimgs(nVols), T2imgs(nVols);
    for (size_t i = 0; i < nVols; i++) {
        inputVector->SetBlockStart(i * multiecho.size());

        apply->SetAlgorithm(algo);
        apply->SetInput(0, inputVector->GetOutput());
        apply->Update();

        PDimgs.at(i) = apply->GetOutput(0);
        T2imgs.at(i) = apply->GetOutput(1);
        PDimgs.at(i)->DisconnectPipeline();
        T2imgs.at(i)->DisconnectPipeline();

        PDoutput->SetInput(i, PDimgs.at(i));
        T2output->SetInput(i, T2imgs.at(i));
    }
    QI_LOG(verbose, "Writing output" );
    PDoutput->UpdateLargestPossibleRegion();
    T2output->UpdateLargestPossibleRegion();
    std::string outPrefix = outarg.Get() + "ME_";
    QI::WriteImage(PDoutput->GetOutput(), outPrefix + "PD" + QI::OutExt());
    QI::WriteImage(T2output->GetOutput(), outPrefix + "T2" + QI::OutExt());
    //QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);

    return EXIT_SUCCESS;
}

