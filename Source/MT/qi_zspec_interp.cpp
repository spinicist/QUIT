/*
 *  qi_cestasym.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>

#include "Util.h"
#include "Args.h"
#include "ImageIO.h"
#include "JSON.h"
#include "Spline.h"
#include "ApplyTypes.h"

class InterpZSpec : public QI::ApplyVectorF::Algorithm {
protected:
    Eigen::ArrayXd m_ifrqs, m_ofrqs;
    TOutput m_zero;
    
public:
    InterpZSpec(const Eigen::ArrayXd &z, const Eigen::ArrayXd &a) :
        m_ifrqs(z), m_ofrqs(a)
    {
        m_zero = TOutput(m_ofrqs.rows()); m_zero.Fill(0.);
    }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_ifrqs.rows(); }
    size_t outputSize() const override {
        return m_ofrqs.rows();
    }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def{{0.0f}};
        return def;
    }
    TOutput zero() const override {
        return m_zero;
    }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"asymmetry"};
        return _names;
    }
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TOutput & /* Unused */,
               TInput & /* Unused */, TIterations & /* Unused */) const override
    {
        const Eigen::Map<const Eigen::ArrayXf> zdata(inputs[0].GetDataPointer(), m_ifrqs.rows());
        QI::SplineInterpolator zspec(m_ifrqs, zdata.cast<double>());
        double f0 = consts.at(0);
        for (int f = 0; f < m_ofrqs.rows(); f++) {
            const double frq = m_ofrqs[f] + f0;
            outputs.at(0)[f] = zspec(frq);
        }
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Interpolates Z-spectrums using Cubic Splines\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPUT", "Change ouput filename (default is input_interp)", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> f0(parser, "OFF RESONANCE", "Specify off-resonance frequency (units must match input)", {'f', "f0"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI_LOG(verbose, "Opening file: " << QI::CheckPos(input_path));
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(input_path));

    QI_LOG(verbose, "Reading input frequences");
    rapidjson::Document json = QI::ReadJSON(std::cin);
    auto i_frqs = QI::ArrayFromJSON(json["input_freqs"]);
    auto o_frqs = QI::ArrayFromJSON(json["output_freqs"]);
    std::shared_ptr<InterpZSpec> algo = std::make_shared<InterpZSpec>(i_frqs, o_frqs);
    auto apply = QI::ApplyVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetInput(0, data);
    if (mask) {
        QI_LOG(verbose, "Setting mask image: " << mask.Get());
        apply->SetMask(QI::ReadImage(mask.Get()));
    }
    if (f0) {
        QI_LOG(verbose, "Setting f0 image: " << f0.Get());
        apply->SetConst(0, QI::ReadImage(f0.Get()));
    }
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(subregion.Get()));
    }
    QI_LOG(verbose, "Processing");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    QI_LOG(verbose, "Elapsed time was " << apply->GetTotalTime() << "s\n");
    std::string outname = outarg ? outarg.Get() : QI::StripExt(input_path.Get()) + "_interp" + QI::OutExt();
    QI_LOG(verbose, "Writing output: " << outname);
    QI::WriteVectorImage(apply->GetOutput(0), outname);
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
