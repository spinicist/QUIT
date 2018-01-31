/*
 *  qi_cestlorentz.cpp
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
#include "ceres/ceres.h"

Eigen::ArrayXd Lorentzian(const double f0, const double fwhm, const double A, const Eigen::ArrayXd &f) {
    Eigen::ArrayXd x = (f0 - f) / (fwhm/2);
    Eigen::ArrayXd s = A / (1. + x.square());
    return s;
}

class ZCost {
private:
    Eigen::ArrayXd m_frqs, m_zspec;
public:

    ZCost(const Eigen::ArrayXd &f, const Eigen::ArrayXd &z) :
          m_frqs(f), m_zspec(z)
    {}

    bool operator() (double const* const* p, double* resids) const {
        const Eigen::ArrayXd L = Lorentzian(p[0][0], p[0][1], p[0][2], m_frqs);
        const Eigen::ArrayXd invL = p[0][3] * (1. - L);
        Eigen::Map<Eigen::ArrayXd> r(resids, m_frqs.size());
        r = (invL - m_zspec);
        return true;
    }
};

class LorentzFit : public QI::ApplyF::Algorithm {
protected:
    Eigen::ArrayXd m_zfrqs;

public:
    LorentzFit(const Eigen::ArrayXf &zf) : m_zfrqs(zf.cast<double>()) { }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 4; }
    size_t dataSize() const override { return m_zfrqs.rows(); }
    size_t outputSize() const override { return 1; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(0, 1.0f);
        return def;
    }
    TOutput zero() const override { return 0; }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"f0", "w", "sat", "PD"};
        return _names;
    }
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TOutput &residual,
               TInput &resids, TIters &its) const override
    {
        const Eigen::Map<const Eigen::ArrayXf> z_spec(inputs[0].GetDataPointer(), m_zfrqs.size());

        // Find closest indices to -2/+2 PPM and only fit Lorentzian between them
        Eigen::ArrayXf::Index indP2, indM2;
        (m_zfrqs + 2.0).abs().minCoeff(&indM2);
        (m_zfrqs - 2.0).abs().minCoeff(&indP2);
        if (indM2 > indP2)
            std::swap(indM2, indP2);
        Eigen::ArrayXf::Index sz = indP2 - indM2;
        const double scale = z_spec.segment(indM2,sz).maxCoeff();
        auto *cost = new ceres::DynamicNumericDiffCostFunction<ZCost>(new ZCost(m_zfrqs.segment(indM2,sz).cast<double>(), z_spec.segment(indM2,sz).cast<double>() / scale));
        cost->AddParameterBlock(4);
        cost->SetNumResiduals(sz);
        Eigen::Array4d p{0.0, 2.0, 0.9, 2.0};
        ceres::Problem problem;
        problem.AddResidualBlock(cost, NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, -2.0);
        problem.SetParameterUpperBound(p.data(), 0, 2.0);
        problem.SetParameterLowerBound(p.data(), 1, 0.001);
        problem.SetParameterUpperBound(p.data(), 1, 100.0);
        problem.SetParameterLowerBound(p.data(), 2, 0.1);
        problem.SetParameterUpperBound(p.data(), 2, 1.0);
        problem.SetParameterLowerBound(p.data(), 3, 0.1);
        problem.SetParameterUpperBound(p.data(), 3, 10.0);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations = 50;
        options.function_tolerance = 1e-5;
        options.gradient_tolerance = 1e-6;
        options.parameter_tolerance = 1e-4;
        ceres::Solve(options, &problem, &summary);
        outputs.at(0) = p[0];
        outputs.at(1) = p[1];
        outputs.at(2) = p[2];
        outputs.at(3) = p[3] * scale;
        residual = summary.final_cost;
        return true;
    }
};


int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Simple Lorentzian fitting.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv);

    if (verbose) std::cout << "Opening file: " << QI::CheckPos(input_path) << std::endl;
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(input_path));

    if (verbose) std::cout << "Enter Z-Spectrum Frequencies: " << std::endl;
    Eigen::ArrayXf z_frqs; QI::ReadArray(std::cin, z_frqs);
    std::shared_ptr<LorentzFit> algo = std::make_shared<LorentzFit>(z_frqs);
    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetInput(0, data);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
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
    std::string outPrefix = outarg.Get() + "LTZ_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        QI::WriteImage(apply->GetOutput(i), outPrefix + algo->names().at(i) + QI::OutExt());
    }
    QI::WriteImage(apply->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
