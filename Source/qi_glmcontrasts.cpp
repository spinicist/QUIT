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

#include "itkTileImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "Filters/MeanImageFilter.h"

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Option.h"

class ContrastsAlgorithm : public QI::ApplyF::Algorithm {
protected:
    Eigen::MatrixXd m_design;
    Eigen::MatrixXd m_contrasts;
    bool m_scale;
public:
    ContrastsAlgorithm(const Eigen::MatrixXd &d, const Eigen::MatrixXd &c, const bool s = false) :
        m_design(d), m_contrasts(c), m_scale(s)
    {}

    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return m_contrasts.rows(); }
    size_t dataSize() const override { return m_design.rows(); }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def;
        return def;
    }
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        Eigen::Map<const Eigen::VectorXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::VectorXd c = m_contrasts * m_design.transpose() * indata.cast<double>() / m_design.rows();
        if (m_scale) {
            c /= indata.mean();
        }
        for (int i = 0; i < m_contrasts.rows(); i++) {
            outputs[i] = c[i];
        }
        residual = 0;
        its = 0;
    }
};


const std::string usage{
"Usage is: qi_glmdiffs input design_matrix contrasts\n\
\n\
A utility for calculating group means, differences etc.\n\
Input is the output file, design matrix and contrasts from qi_glmsetup.\n\
One output file will be generated for each contrast.\n"};

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts(usage);
    QI::Switch fractional('F',"frac","Output contrasts as fraction of grand mean", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Add a prefix to output filenames", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 3) {
        std::cerr << opts << std::endl;
        std::cerr << "Must specify input, design_matrix and contrasts" << std::endl;
        return EXIT_FAILURE;
    }

    if (*verbose) std::cout << "Reading input file" << std::endl;
    QI::VectorVolumeF::Pointer merged = QI::ReadVectorImage<float>(nonopts.at(0));
    if (*verbose) std::cout << "Reading design matrix" << std::endl;
    std::ifstream design_matrix_file(nonopts.at(1));
    Eigen::ArrayXXd design_matrix;
    QI::ReadArray(design_matrix_file, design_matrix);
    if (*verbose) std::cout << "Reading contrasts file" << std::endl;
    std::ifstream contrasts_file(nonopts.at(2));
    Eigen::ArrayXXd contrasts;
    QI::ReadArray(contrasts_file, contrasts);

    if (design_matrix.rows() != merged->GetNumberOfComponentsPerPixel()) {
        std::cerr << "Number of rows in design matrix (" << design_matrix.rows()
                  << ") does not match number of volumes in image (" << merged->GetNumberOfComponentsPerPixel() << ")" << std::endl;
        return EXIT_FAILURE;
    }

    auto con_algo = std::make_shared<ContrastsAlgorithm>(design_matrix, contrasts, *fractional);
    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(con_algo);
    apply->SetInput(0, merged);
    apply->SetMask(*mask);
    if (*verbose) std::cout << "Calculating contrasts" << std::endl;
    apply->Update();
    for (int c = 0; c < contrasts.rows(); c++) {
        if (*verbose) std::cout << "Writing contrast " << c << std::endl;
        QI::WriteImage(apply->GetOutput(c), *outPrefix + "con" + std::to_string(c) + ".nii");
    }
}

