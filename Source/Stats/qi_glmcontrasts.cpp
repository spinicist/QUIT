/*
 *  qi_glmdiffs.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <Eigen/Dense>

#include "itkTileImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "MeanImageFilter.h"

#include "ApplyTypes.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

class ContrastsAlgorithm : public QI::ApplyF::Algorithm {
protected:
    Eigen::MatrixXd m_mat;
    bool m_scale;
public:
    ContrastsAlgorithm(const Eigen::MatrixXd &d, const Eigen::MatrixXd &c, const bool s = false) :
        m_scale(s)
    {
        m_mat = c * ((d.transpose() * d).inverse()) * d.transpose();
    }

    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return m_mat.rows(); }
    size_t dataSize() const override { return m_mat.cols(); }
    float zero() const override { return 0.f; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def;
        return def;
    }
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::VectorXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::VectorXd c = m_mat * indata.cast<double>();
        if (m_scale) {
            c /= indata.mean();
        }
        for (int i = 0; i < m_mat.rows(); i++) {
            outputs[i] = c[i];
        }
        residual = 0;
        its = 0;
        return true;
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("A utility for calculating group means, differences etc.\n"
                                "One output file will be generated for each contrast.\n"
                                "\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "IMAGE", "The combined image file from qi_glmsetup");
    args::Positional<std::string> design_path(parser, "DESIGN", "GLM Design matrix from qi_glmsetup");
    args::Positional<std::string> contrasts_path(parser, "CONTRASTS", "Contrasts matrix from qi_glmsetup");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag fraction(parser, "FRACTION", "Output contrasts as fraction of grand mean", {'F',"frac"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    QI::ParseArgs(parser, argc, argv);

    if (verbose) std::cout << "Reading input file " << QI::CheckPos(input_path) << std::endl;
    QI::VectorVolumeF::Pointer merged = QI::ReadVectorImage<float>(QI::CheckPos(input_path));
    if (verbose) std::cout << "Reading design matrix" << QI::CheckPos(design_path) << std::endl;
    Eigen::ArrayXXd design_matrix = QI::ReadArrayFile(QI::CheckPos(design_path));
    if (verbose) std::cout << "Reading contrasts file" << QI::CheckPos(contrasts_path) << std::endl;
    Eigen::ArrayXXd contrasts     = QI::ReadArrayFile(QI::CheckPos(contrasts_path));
    if (design_matrix.rows() != merged->GetNumberOfComponentsPerPixel()) {
        std::cerr << "Number of rows in design matrix (" << design_matrix.rows()
                  << ") does not match number of volumes in image (" << merged->GetNumberOfComponentsPerPixel() << ")" << std::endl;
        return EXIT_FAILURE;
    }
    if (design_matrix.cols() != contrasts.cols()) {
        std::cerr << "Number of columns in design matrix (" << design_matrix.cols()
                  << ") does not match contrasts (" << contrasts.cols() << ")" << std::endl;
        return EXIT_FAILURE;
    }

    auto con_algo = std::make_shared<ContrastsAlgorithm>(design_matrix, contrasts, fraction);
    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(con_algo);
    apply->SetInput(0, merged);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (verbose) std::cout << "Calculating contrasts" << std::endl;
    apply->Update();
    for (int c = 0; c < contrasts.rows(); c++) {
        if (verbose) std::cout << "Writing contrast " << (c + 1) << std::endl;
        QI::WriteImage(apply->GetOutput(c), outarg.Get() + "con" + std::to_string(c + 1) + QI::OutExt());
    }
}

