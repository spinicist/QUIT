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
#include "IO.h"
#include "Spline.h"

class MTAsym : public QI::ApplyVectorF::Algorithm {
protected:
    Eigen::ArrayXd m_zfrqs, m_afrqs;
    TOutput m_zero;
    
public:
    MTAsym(const Eigen::ArrayXf &zf, const Eigen::ArrayXf &af) :
        m_zfrqs(zf.cast<double>()), m_afrqs(af.cast<double>())
    {
        m_zero = TOutput(m_afrqs.rows()); m_zero.Fill(0.);
    }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_zfrqs.rows(); }
    size_t outputSize() const override {
        return m_afrqs.rows();
    }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 0.0f);
        return def;
    }
    TOutput zero() const override {
        return m_zero;
    }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"asymmetry"};
        return _names;
    }
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TOutput &residual,
               TInput &resids, TIters &its) const override
    {
        const Eigen::Map<const Eigen::ArrayXf> zdata(inputs[0].GetDataPointer(), m_zfrqs.rows());
        QI::SplineInterpolator zspec(m_zfrqs, zdata.cast<double>());
        double ref = zdata[0];
        double f0 = consts.at(0);
        for (int f = 0; f < m_afrqs.rows(); f++) {
            const double pfrq = m_afrqs[f] - f0;
            const double nfrq = f0 - m_afrqs[f];
            const double pos = zspec(pfrq);
            const double neg = zspec(nfrq);
            outputs.at(0)[f] = ((pos - neg)/ref);
        }
        return true;
    }
};


int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Simple MT-asymmetry calculator.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> f0(parser, "OFF RESONANCE", "Specify off-resonance frequency", {'f', "f0"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;

    if (verbose) std::cout << "Opening file: " << QI::CheckPos(input_path) << std::endl;
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(input_path));

    if (prompt) std::cout << "Enter Z-Spectrum Frequencies: " << std::endl;
    Eigen::ArrayXf z_frqs; QI::ReadArray(std::cin, z_frqs);
    if (prompt) std::cout << "Enter Asymmetry Frequencies: " << std::endl;
    Eigen::ArrayXf a_frqs; QI::ReadArray(std::cin, a_frqs); // Asymmetry output
    std::shared_ptr<MTAsym> algo = std::make_shared<MTAsym>(z_frqs, a_frqs);
    auto apply = QI::ApplyVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetInput(0, data);
    if (mask) {
        if (verbose) std::cout << "Setting mask image: " << mask.Get() << std::endl;
        apply->SetMask(QI::ReadImage(mask.Get()));
    }
    if (f0) {
        if (verbose) std::cout << "Setting f0 image: " << f0.Get() << std::endl;
        apply->SetConst(0, QI::ReadImage(f0.Get()));
    }
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(subregion.Get()));
    }
    if (verbose) {
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing output." << std::endl;
    }
    std::string outPrefix = outarg.Get() + "CEST_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        QI::WriteVectorImage(apply->GetOutput(i), outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
