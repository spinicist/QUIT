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
#include "IO.h"
#include "Spline.h"
#include "ApplyTypes.h"

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
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TOutput & /* Unused */,
               TInput & /* Unused */, TIterations & /* Unused */) const override
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
        return std::make_tuple(true, "");
    }
};


int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Simple MT-asymmetry calculator.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> f0(parser, "OFF RESONANCE", "Specify off-resonance frequency", {'f', "f0"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);

    QI_LOG(verbose, "Opening file: " << QI::CheckPos(input_path));
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(input_path));

    QI_LOG(verbose, "Enter Z-Spectrum Frequencies: " );
    Eigen::ArrayXf z_frqs; QI::ReadArray(std::cin, z_frqs);
    QI_LOG(verbose, "Enter Asymmetry Frequencies: " );
    Eigen::ArrayXf a_frqs; QI::ReadArray(std::cin, a_frqs); // Asymmetry output
    std::shared_ptr<MTAsym> algo = std::make_shared<MTAsym>(z_frqs, a_frqs);
    auto apply = QI::ApplyVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
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
    QI_LOG(verbose, "Elapsed time was " << apply->GetTotalTime() << "s" <<
                    "Writing output.");
    std::string outPrefix = outarg.Get() + "MT_";
    for (size_t i = 0; i < algo->numOutputs(); i++) {
        QI::WriteVectorImage(apply->GetOutput(i), outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
