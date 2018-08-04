/*
 *  qi_zshim.cpp
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
#include "ImageIO.h"
#include "Args.h"
#include "ApplyTypes.h"

class ZShimSoS : public QI::ApplyVectorF::Algorithm {
protected:
    const int m_inputsize, m_zshims;
    bool m_debug = false;
    TOutput m_zero;
public:
    ZShimSoS(const int insize, const int zshims, bool debug) :
        m_inputsize(insize), m_zshims(zshims), m_debug(debug)
    {
        m_zero = TOutput(this->outputSize());
        m_zero.Fill(0.);
    }

    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 1; }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"zshim"};
        return _names;
    }
    size_t dataSize() const override { return m_inputsize; }
    size_t outputSize() const override { return m_inputsize / m_zshims; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(0);
        return def;
    }
    TOutput zero() const override { return m_zero; }
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> &/* Unused */,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TOutput &/* Unused */,
               TInput &/* Unused */, TIterations &/* Unused */) const override
    {
        Eigen::Map<const Eigen::MatrixXf> input(inputs[0].GetDataPointer(), m_zshims, m_inputsize / m_zshims);
        Eigen::VectorXf output = input.colwise().norm();
        for (size_t i = 0; i < this->outputSize(); i++) {
            outputs[0][i] = output[i];
        }
        return std::make_tuple(true, "");
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Combines Z-Shimmed data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ZSHIM_FILE", "Input Z-Shimmed file");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Output progress messages", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debug messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<int> zshims(parser, "ZSHIMS", "Number of Z-Shims (default 8)", {'z', "zshims"}, 8);
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);
    QI_LOG(verbose, "Reading Z-Shim data from: " << QI::CheckPos(input_path));
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path));
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    auto algo = std::make_shared<ZShimSoS>(input->GetNumberOfComponentsPerPixel(), zshims.Get(), debug);
    auto apply = QI::ApplyVectorF::New();
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(false);
    QI_LOG(verbose, "Using " << threads.Get() << " threads" );
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    apply->SetInput(0, input);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    }
    QI_LOG(verbose, "Processing");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    QI_LOG(verbose, "Elapsed time was " << apply->GetTotalTime() << "s");
    const std::string fname = outPrefix + "_" + algo->names()[0] + QI::OutExt();
    QI_LOG(verbose, "Writing file: " << fname);
    QI::WriteVectorImage(apply->GetOutput(0), fname);
    return EXIT_SUCCESS;
}
