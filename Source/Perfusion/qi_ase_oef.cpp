/*
 *  qi_ase_oef.cpp
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

#include "itkDerivativeImageFilter.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "ApplyTypes.h"
#include "MultiEchoSequence.h"
#include "SequenceCereal.h"

class ASEAlgo : public QI::ApplyF::Algorithm {
protected:
    const std::shared_ptr<QI::MultiEchoBase> m_sequence;
    const int m_inputsize;
    const double m_B0;
    double m_Tc;
    const QI::VolumeF::SpacingType m_voxsize;
    Eigen::ArrayXd m_TE; // Local copy excluding values below Tc
    Eigen::Index m_TE0;
    // Constants for calculations
    const double kappa = 0.03; // Conversion factor
    const double gamma = 42.577e6; // Gyromagnetic Ratio
    const double delta_X0 = 0.264e-6; // Difference in susceptibility of oxy and fully de-oxy blood
    const double Hb = 0.34 / kappa; // Hct = 0.34;

public:
    ASEAlgo(const std::shared_ptr<QI::MultiEchoBase> &seq, const int inputsize, const double B0, const QI::VolumeF::SpacingType voxsize) :
        m_sequence(seq), m_inputsize(inputsize), m_B0(B0), m_voxsize(voxsize)
    {
        // Nic Blockley uses Tc = 15 ms for 3T, scale for other field-strengths
        m_Tc = 0.015 / (B0 / 3);
        const int above_Tc_count = (m_sequence->TE.abs() > m_Tc).count();
        m_TE = Eigen::ArrayXd(above_Tc_count);
        Eigen::Index ind = 0;
        m_TE0 = m_sequence->size();
        for (Eigen::Index i = 0; i < m_sequence->size(); i++) {
            if (m_sequence->TE[i] == 0) {
                m_TE0 = i;
            }
            if (std::abs(m_sequence->TE[i]) > m_Tc) {
                m_TE[ind] = std::abs(m_sequence->TE[i]);
                ind++;
            }
        }
        if (m_TE0 == m_sequence->size()) {
            QI_FAIL("Did not find a zero echo-time in input");
        }
        // QI_DB( B0 );
        // QI_DB( m_Tc );
        // QI_DB( m_TE0 );
        // QI_DB( above_Tc_count );
        // QI_DBVEC( m_sequence->TE );
        // QI_DBVEC( m_TE );
    }

    size_t numInputs() const override  { return 1; }
    size_t numConsts() const override  { return 3; }
    size_t numOutputs() const override { return 4; }
    size_t dataSize() const override   { return m_inputsize; }
    size_t outputSize() const override { return 1; }
    TOutput zero() const override {
        TOutput z{0};
        return z;
    }

    std::vector<float> defaultConsts() const override {
        std::vector<float> def(3, 0.0); // No field gradients
        return def;
    }

    const std::vector<std::string> &names() const {
        static std::vector<std::string> _names = {"R2prime", "DBV", "OEF", "dHb"};
        return _names;
    }

    double sinc(const double x) const {
        static double const taylor_0_bound = std::numeric_limits<double>::epsilon();
        static double const taylor_2_bound = sqrt(taylor_0_bound);
        static double const taylor_n_bound = sqrt(taylor_2_bound);

        if (std::abs(x) >= taylor_n_bound) {
            return(sin(x)/x);
        } else {
            // approximation by taylor series in x at 0 up to order 0
            double result = 1;

            if (abs(x) >= taylor_0_bound) {
                double x2 = x*x;
                // approximation by taylor series in x at 0 up to order 2
                result -= x2/6;

                if (abs(x) >= taylor_2_bound) {
                    // approximation by taylor series in x at 0 up to order 4
                    result += (x2*x2)/120;
                }
            }
            return result;
        }
    }
    
    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex & /* Unused */,
               std::vector<TOutput> &outputs, TOutput &residual,
               TInput &resids, TIterations &its) const override
    {
        const Eigen::Map<const Eigen::ArrayXXf> input(inputs[0].GetDataPointer(), m_sequence->size(), inputs[0].Size() / m_sequence->size());
        Eigen::ArrayXd all_data = input.cast<double>().rowwise().mean();
        Eigen::ArrayXd data(m_TE.rows());
        Eigen::Index ind = 0;
        for (Eigen::Index i = 0; i < all_data.rows(); i++) {
            if (std::abs(m_sequence->TE[i]) > m_Tc) {
                double F = 1.0;
                for (auto d = 0; d < 3; d++) {
                    const double grad = consts[d];
                    const double x = grad * 2 * M_PI * m_voxsize[d] * m_sequence->TE[i] / 2;
                    F *= std::abs(sinc(x));
                }
                data[ind] = all_data[i] / F;
                ind++;
            }
        }

        Eigen::MatrixXd X(m_TE.rows(), 2);
        X.col(0) = m_TE;
        X.col(1).setOnes();
        Eigen::VectorXd Y = data.log();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        const double R2prime = -b[0];
        const double logS0_linear = b[1];
        const double DBV = logS0_linear - log(data[m_TE0]);
        const double dHb = 3*R2prime / (DBV * 4 * gamma * M_PI * delta_X0 * kappa * m_B0);
        const double OEF = dHb / Hb;

        outputs[0] = R2prime;
        outputs[1] = DBV*100;
        outputs[2] = OEF*100;
        outputs[3] = dHb;
        residual = 0;
        resids.Fill(0.);
        its = 0;
        return std::make_tuple(true, "");
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates the OEF from ASE data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ASE_FILE", "Input ASE file");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<double> B0(parser, "B0", "Field-strength (Tesla), default 3", {'B', "B0"}, 3.0);
    args::ValueFlag<std::string> gradx(parser, "GRADX", "Gradient of field-map in x-direction for MFG correction", {'x', "gradx"});
    args::ValueFlag<std::string> grady(parser, "GRADY", "Gradient of field-map in y-direction for MFG correction", {'y', "grady"});
    args::ValueFlag<std::string> gradz(parser, "GRADZ", "Gradient of field-map in z-direction for MFG correction", {'z', "gradz"});
    args::ValueFlag<double> slice_arg(parser, "SLICE THICKNESS", "Slice-thickness for MFG calculation (useful if there was a slice gap)", {'s', "slice"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);
    if (verbose) std::cout << "Reading ASE data from: " << QI::CheckPos(input_path) << std::endl;
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path));
    std::shared_ptr<QI::MultiEchoBase> sequence;
    cereal::JSONInputArchive ar(std::cin);
    load(ar, sequence);
    QI::VolumeF::SpacingType vox_size = input->GetSpacing();
    if (slice_arg) {
        vox_size[2]  = slice_arg.Get();
    }
    std::shared_ptr<ASEAlgo> algo = std::make_shared<ASEAlgo>(sequence, input->GetNumberOfComponentsPerPixel(), B0.Get(), vox_size);
    auto apply = QI::ApplyF::New();
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(false);
    if (verbose) std::cout << "Using " << threads.Get() << " threads" << std::endl;
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    apply->SetInput(0, input);
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (gradx) apply->SetConst(0, QI::ReadImage(gradx.Get()));
    if (grady) apply->SetConst(1, QI::ReadImage(grady.Get()));
    if (gradz) apply->SetConst(2, QI::ReadImage(gradz.Get()));
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    }
    if (verbose) {
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
    }
    
    for (size_t i = 0; i < algo->numOutputs(); i++) {
        const std::string fname = outPrefix + "_" + algo->names()[i] + QI::OutExt();
        std::cout << "Writing file: " << fname << std::endl;
        QI::WriteImage(apply->GetOutput(i), fname);
    }
    return EXIT_SUCCESS;
}
