/*
 *  qdespot1.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include <Eigen/Dense>
#include "ceres/ceres.h"

#include "ApplyTypes.h"
#include "SSFPSequence.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

struct PLANET : public QI::ApplyVectorF::Algorithm {
    const QI::SSFPGSSequence &m_seq;

    size_t numInputs() const override { return 3; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override { return m_seq.FA.rows() * 3; }
    size_t outputSize() const override { return m_seq.FA.rows(); }
    TOutput zero() const override {
        TOutput zero(m_seq.FA.rows());
        zero.Fill(0.);
        return zero;
    }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 1.0); // B1
        return def;
    }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"T1", "T2", "PD"};
        return _names;
    }
    PLANET(const QI::SSFPGSSequence &s) : m_seq(s) {}

    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &/*Unused*/,
               std::vector<TOutput> &outputs, TOutput &/*Unused*/,
               TInput &/*Unused*/, TIterations &/*Unused*/) const override
    {
        Eigen::Map<const Eigen::ArrayXf> G(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::Map<const Eigen::ArrayXf> a(inputs[1].GetDataPointer(), inputs[0].Size());
        Eigen::Map<const Eigen::ArrayXf> b(inputs[2].GetDataPointer(), inputs[0].Size());
        const float b1 = consts[0];
        const Eigen::ArrayXf cosa = cos(b1 * m_seq.FA.cast<float>());
        const Eigen::ArrayXf sina = sin(b1 * m_seq.FA.cast<float>());
        const Eigen::ArrayXf T1 = -m_seq.TR / log((a*(1. + cosa - a*b*cosa) - b)/(a*(1. + cosa - a*b) - b*cosa));
        const Eigen::ArrayXf T2 = -m_seq.TR / log(a);
        const Eigen::ArrayXf E1 = exp(-m_seq.TR / T1);
        const Eigen::ArrayXf E2 = a; // For simplicity copying formulas
        const Eigen::ArrayXf PD = G * (1. - E1*cosa - E2*E2*(E1 - cosa)) / (sqrt(E2)*(1. - E1)*sina);
        for (int i = 0; i < m_seq.FA.rows(); i++) {
            outputs[0][i] = T1[i];
            outputs[1][i] = T2[i];
            outputs[2][i] = PD[i];
        }
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    args::ArgumentParser parser("Calculates T1&T2 from SSFP Ellipse Parameters.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> G_filename(parser, "G", "Ellipse parameter G");
    args::Positional<std::string> a_filename(parser, "a", "Ellipse parameter a");
    args::Positional<std::string> b_filename(parser, "b", "Ellipse parameter b");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> out_prefix(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threads.Get());
    QI_LOG(verbose, "Opening G: " << QI::CheckPos(G_filename));
    auto G = QI::ReadVectorImage(QI::CheckPos(G_filename));
    QI_LOG(verbose, "Opening a: " << QI::CheckPos(a_filename));
    auto a = QI::ReadVectorImage(QI::CheckPos(a_filename));
    QI_LOG(verbose, "Opening b: " << QI::CheckPos(b_filename));
    auto b = QI::ReadVectorImage(QI::CheckPos(b_filename));
    rapidjson::Document json = QI::ReadJSON(std::cin);
    QI::SSFPGSSequence seq(json["SSFPGS"]);
    auto algo = std::make_shared<PLANET>(seq);
    auto apply = QI::ApplyVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetInput(0, G);
    apply->SetInput(1, a);
    apply->SetInput(2, b);
    if (B1) apply->SetConst(0, QI::ReadImage(B1.Get()));
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
    std::string outPrefix = out_prefix.Get() + "PLANET_";
    for (size_t i = 0; i < algo->numOutputs(); i++) {
        QI_LOG(verbose, "Writing output: " << outPrefix + algo->names().at(i) + QI::OutExt());
        QI::WriteVectorImage(apply->GetOutput(i), outPrefix + algo->names().at(i) + QI::OutExt());
    }
    QI_LOG(verbose, "Finished." );
    return EXIT_SUCCESS;
}
