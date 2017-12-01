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
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Splines>

#include "QI/Util.h"
#include "QI/Args.h"
#include "QI/IO.h"

class MTAsym : public QI::ApplyVectorF::Algorithm {
protected:
    Eigen::ArrayXf m_zfrqs, m_afrqs;
    TOutput m_zero;
public:
    MTAsym(const Eigen::ArrayXf &zf, const Eigen::ArrayXf &af) :
        m_zfrqs(zf), m_afrqs(af)
    {
        m_zero = TOutput(m_afrqs.rows()); m_zero.Fill(0.);
    }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_zfrqs.rows(); }
    size_t outputSize() const override {
        return m_afrqs.rows();
    }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(0, 1.0f);
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
        typedef Eigen::Spline<float, 1> TSpline;
        typedef Eigen::SplineFitting<TSpline> TFit;
        const float maxfrq = m_zfrqs.maxCoeff();
        const float minfrq = m_zfrqs.minCoeff();
        const float w = maxfrq - minfrq;
        const Eigen::ArrayXf scaledfrqs = (m_zfrqs - minfrq) / w;
        Eigen::DenseIndex degree = std::min<int>(m_zfrqs.rows() - 1, 3);
        const Eigen::Map<const Eigen::ArrayXf> zdata(inputs[0].GetDataPointer(), m_zfrqs.rows());
        // std::cout << "zfrqs " << m_zfrqs.transpose() << std::endl;
        // std::cout << "sfrqs " << scaledfrqs.transpose() << std::endl;
        // std::cout << "afrqs " << m_afrqs.transpose() << std::endl; 
        // std::cout << "zdata " << zdata.transpose() << std::endl;
        // std::cout << "degree " << degree << std::endl;
        // std::cout << "zdata " << zdata.transpose() << std::endl;
        // std::cout << "degree " << degree << std::endl;
        TSpline zspec;
        if (scaledfrqs[0] > 0) {
            // std::cout << "Reverse" << std::endl;
            zspec = TFit::Interpolate(zdata.reverse().transpose(), degree, scaledfrqs.reverse());
        } else {
            // std::cout << "Forward" << std::endl;
            zspec = TFit::Interpolate(zdata.transpose(), degree, scaledfrqs);
        }
        float ref = zdata[0];
        for (int f = 0; f < m_afrqs.rows(); f++) {
            const float pfrq =  (m_afrqs.coeffRef(f) - minfrq)/w;
            const float nfrq = (-m_afrqs.coeffRef(f) - minfrq)/w;
            const float pos = zspec(pfrq)[0];
            const float neg = zspec(nfrq)[0];
            // std::cout << pos << " " << neg << " " << ref << " " << ((pos - neg)/ref) << std::endl;
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
    apply->SetMask(QI::ReadImage(mask.Get()));
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
