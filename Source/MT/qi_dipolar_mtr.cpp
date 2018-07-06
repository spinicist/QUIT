/*
 *  qi_dipolar_mtr.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>

#include "ApplyTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"

class DMTR : public QI::ApplyF::Algorithm {
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 4; }
    size_t dataSize() const override { return 5; }
    float zero() const override { return 0.f; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def;
        return def;
    }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"mtr", "emtr", "dmtr", "mta"};
        return _names;
    }

    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> & /* Unused */,
               const TIndex & /* Unused */,
               std::vector<TOutput> &outputs, TConst & /* Unused */,
               TInput & /* Unused */, TIterations & /* Unused */) const override
    {
        outputs[0] = 100. * (1. - (inputs[0][3] + inputs[0][4]) / (2. * inputs[0][2]));
        outputs[1] = 100. * (1. - (inputs[0][0] + inputs[0][1]) / (2. * inputs[0][2]));
        outputs[2] = outputs[1] - outputs[0];
        outputs[3] = 100. * (inputs[0][3] - inputs[0][4]) / inputs[0][2];
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    args::ArgumentParser parser("Calculates MTR & dipolar MTR.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_file(parser, "MT_FILE", "Input image. Must have 5 volumes (Sat+, Sat+-, Unsat, Sat-+, Sat-)");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> out_prefix(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);

    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threads.Get());
    QI_LOG(verbose, "Opening MT file " << QI::CheckPos(input_file));
    auto volumes = QI::ReadVectorImage(QI::CheckPos(input_file));
    auto algo = std::make_shared<DMTR>();
    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetInput(0, volumes);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(subregion.Get()));
    }
    QI_LOG(verbose, "Processing");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    QI_LOG(verbose, "Elapsed time was " << apply->GetTotalTime() << "s");
    std::string outPrefix = out_prefix.Get() + "DMT_";
    for (size_t i = 0; i < algo->numOutputs(); i++) {
        QI_LOG(verbose, "Writing output: " << outPrefix + algo->names().at(i) + QI::OutExt());
        QI::WriteImage(apply->GetOutput(i), outPrefix + algo->names().at(i) + QI::OutExt());
    }
    QI_LOG(verbose, "Finished." );
    return EXIT_SUCCESS;
}

