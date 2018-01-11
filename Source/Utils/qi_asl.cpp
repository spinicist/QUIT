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

#include <iostream>
#include <Eigen/Dense>

#include "Util.h"
#include "IO.h"
#include "Args.h"
#include "Models.h"
#include "Sequences.h"
#include "Types.h"

class CASLSequence {
public:
    double label_time, post_label_delay, TR;

    CASLSequence(std::istream& istr, const bool prompt) {
        if (prompt) std::cout << "Enter TR: " << std::flush;
        QI::Read(istr, this->TR);
        if (prompt) std::cout << "Enter label time: " << std::flush;
        QI::Read(istr, this->label_time);
        if (prompt) std::cout << "Enter post-label delay: " << std::flush;
        QI::Read(istr, this->post_label_delay);
    }
    void write(std::ostream &os) const {
        os << "CASL" << std::endl;
        os << "TR: " << this->TR
           << "\tLabel Time: " << this->label_time
           << "\tPost-Label Delay: " << this->post_label_delay << std::endl;
    }
};

class CASLAlgo : public QI::ApplyVectorF::Algorithm {
protected:
    const CASLSequence m_CASL;
    const double m_T1, m_alpha, m_lambda;
    const int m_inputsize, m_series_size;
    const bool m_average_timeseries;
public:
    CASLAlgo(const CASLSequence& casl,
             const double T1, const double alpha, const double lambda,
             const int inputsize, const bool average) :
        m_CASL(casl), m_T1(T1), m_alpha(alpha), m_lambda(lambda),
        m_inputsize(inputsize), m_series_size(inputsize/2),
        m_average_timeseries(average)
    {
    }

    size_t numInputs() const override  { return 1; }
    size_t numConsts() const override  { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override   { return m_inputsize; }
    size_t outputSize() const override {
        if (m_average_timeseries) {
            return 1;
        } else {
            return m_series_size;
        }
    }
    TOutput zero() const override {
        TOutput z;
        if (m_average_timeseries) {
            z.SetSize(1);
        } else {
            z.SetSize(m_series_size);
        }
        z.Fill(0.);
        return z;
    }

    std::vector<float> defaultConsts() const override {
        std::vector<float> def(0, 0);
        return def;
    }

    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TOutput &residual,
               TInput &resids, TIters &its) const override
    {
        const Eigen::Map<const Eigen::ArrayXf, 0, Eigen::InnerStride<>> even(inputs[0].GetDataPointer(), m_series_size, Eigen::InnerStride<>(2));
        const Eigen::Map<const Eigen::ArrayXf, 0, Eigen::InnerStride<>> odd(inputs[0].GetDataPointer() + 1, m_series_size, Eigen::InnerStride<>(2));

        const Eigen::ArrayXd diff = (odd.cast<double>() - even.cast<double>());
        const Eigen::ArrayXd SI_PD = odd.cast<double>();
        const Eigen::ArrayXd CBF = (6000 * m_lambda * diff * exp(m_CASL.post_label_delay / m_T1)) / 
                           (2. * m_alpha * m_T1 * SI_PD * (1. - exp(-m_CASL.label_time / m_T1)));
        // std::cout << "l " << m_lambda << " diff " << diff << " PLD " << m_PLD << " T1 " << m_T1 << " e(PLD) " << exp(m_PLD / m_T1) << std::endl;
        // std::cout << "a " << m_alpha << " PD " << SI_PD << " LD " << m_LD << " (1 - exp()) " << (1. - exp(-m_LD / m_T1)) << std::endl;
        // std::cout << "CBF: " << CBF << std::endl;
        if (m_average_timeseries) {
            outputs[0][0] = CBF.mean();
        } else {
            for (int i = 0; i < m_series_size; i++) {
                outputs[0][i] = CBF[i];
            }
        }
        residual.Fill(0.);
        resids.Fill(0.);
        its = 0;
        return true;
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates CBF from ASL data.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "ASL_FILE", "Input ASL file");
    
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::Flag              average(parser, "AVERAGE", "Average the time-series", {'a', "average"});
    args::ValueFlag<double> T1_blood(parser, "BLOOD T1", "Value of blood T1 to use (seconds), default 2.429", {'t', "T1"}, 2.429);
    args::ValueFlag<double> alpha(parser, "ALPHA", "Labelling efficiency, default 0.9", {'a', "alpha"}, 0.9);
    args::ValueFlag<double> lambda(parser, "LAMBDA", "Blood-brain partition co-efficent, default 0.9 mL/g", {'l', "lambda"}, 0.9);
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;

    if (verbose) std::cout << "Reading ASL data from: " << QI::CheckPos(input_path) << std::endl;
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path));

    CASLSequence sequence(std::cin, prompt);

    auto apply = QI::ApplyVectorF::New();
    if (verbose) {
        sequence.write(std::cout);
    }
    std::shared_ptr<CASLAlgo> algo = std::make_shared<CASLAlgo>(sequence, T1_blood.Get(), alpha.Get(), lambda.Get(),
                                                                input->GetNumberOfComponentsPerPixel(), average);
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(false);
    if (verbose) std::cout << "Using " << threads.Get() << " threads" << std::endl;
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    apply->SetInput(0, input);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(args::get(subregion)));
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
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    QI::WriteVectorImage(apply->GetOutput(0), outPrefix + "CBF" + QI::OutExt());
    return EXIT_SUCCESS;
}
