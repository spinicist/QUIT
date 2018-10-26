/*
 *  qimp2rage.cpp
 *
 *  Created by Tobias Wood on 2015/08/24.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <complex>

#include "itkUnaryGeneratorImageFilter.h"
#include "itkBinaryGeneratorImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkAddImageFilter.h"

// #define QI_DEBUG_BUILD 1
#include "ImageTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "MPRAGESequence.h"
#include "Masking.h"
#include "Spline.h"

inline float MP2Contrast(const std::complex<float> ti1, 
                         const std::complex<float> ti2,
                         const float beta = 0.0f)
{
    const float a1 = abs(ti1);
    const float a2 = abs(ti2);
    return (std::real(std::conj(ti1)*ti2) - beta)/(std::norm(a1) + std::norm(a2) + 2*beta);
}

Eigen::Array2cf One_MP2RAGE(const double &M0, const double &T1, const double &B1, const QI::MP2RAGESequence &s) {
    const double R1 = 1. / T1;
    const Eigen::Array2d R1s = R1 - log(cos(B1 * s.FA))/s.TR;
    const Eigen::Array2d M0s = M0 * (1. - exp(-s.TR*R1)) / (1. - exp(-s.TR*R1s));
    const double tau = s.SegLength * s.TR;

    const Eigen::Array3d B = exp(-s.TD*R1);
    const Eigen::Array3d A = M0*(1. - B);

    const Eigen::Array2d D = exp(-tau*R1s);
    const Eigen::Array2d C = M0s*(1. - D);

    Eigen::Array2d Mm;
    const double eta = 1.0;
    const double denominator = (1 + eta*B[0]*D[0]*B[1]*D[1]*B[2]);
    Mm[0] = (A[0]-eta*B[0]*(A[2]+B[2]*(C[1]+D[1]*(A[1]+B[1]*C[0])))) / denominator;
    Mm[1] = (A[1]+B[1]*(C[0]+D[0]*(A[0]-eta*B[0]*(A[2]+B[2]*C[1])))) / denominator;
    Eigen::Array2cf Me = Eigen::Array2cf::Zero();
    Me.real() = ((M0s + (Mm - M0s)*exp(-s.TR*R1s*s.k0))*sin(B1 * s.FA)).cast<float>();
    QI_DB( M0 )
    QI_DB( T1 )
    QI_DB( B1 )
    QI_DBVEC( Mm )
    QI_DBVEC( Me )
    QI_DB( s.TR )
    QI_DBVEC( s.TD )
    QI_DBVEC( s.FA )
    QI_DB( s.SegLength )
    QI_DB( s.k0 )
    return Me;
}

int main(int argc, char **argv) {
    args::ArgumentParser parser("Calculates T1/B1 maps from MP2/3-RAGE data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT FILE", "Path to complex MP-RAGE data");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> seq_arg(parser, "FILE", "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float> beta_arg(parser, "BETA", "Regularisation factor for robust contrast calculation (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676)", {'b', "beta"}, 0.0);
    args::Flag t1(parser, "T1", "Calculate T1 map via spline look-up", {'t', "t1"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI_LOG(verbose, "Opening input file " << QI::CheckPos(input_path));
    auto inFile = QI::ReadImage<QI::SeriesXF>(QI::CheckPos(input_path));

    auto ti_1 = itk::ExtractImageFilter<QI::SeriesXF, QI::VolumeXF>::New();
    auto ti_2 = itk::ExtractImageFilter<QI::SeriesXF, QI::VolumeXF>::New();
    auto region = inFile->GetLargestPossibleRegion();
    region.GetModifiableSize()[3] = 0;
    ti_1->SetExtractionRegion(region);
    ti_1->SetDirectionCollapseToSubmatrix();
    ti_1->SetInput(inFile);
    region.GetModifiableIndex()[3] = 1;
    ti_2->SetExtractionRegion(region);
    ti_2->SetDirectionCollapseToSubmatrix();
    ti_2->SetInput(inFile);

    QI_LOG(verbose, "Generating MP2 contrasts");
    using BinaryFilter = itk::BinaryGeneratorImageFilter<QI::VolumeXF, QI::VolumeXF, QI::VolumeF>;
    auto MP2Filter = BinaryFilter::New();
    MP2Filter->SetInput1(ti_1->GetOutput());
    MP2Filter->SetInput2(ti_2->GetOutput());
    const float &beta = beta_arg.Get();
    MP2Filter->SetFunctor( [&](const std::complex<float> &p1, const std::complex<float> &p2) { return MP2Contrast(p1, p2, beta); });
    MP2Filter->Update();
    const std::string out_prefix = outarg ? outarg.Get() : QI::StripExt(input_path.Get());
    QI::WriteImage(MP2Filter->GetOutput(), out_prefix + "_MP2" + QI::OutExt(), verbose);

    if (t1) {
        QI_LOG(verbose, "Reading sequence information");
        rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
        QI::MP2RAGESequence mp2rage_sequence(input["MP2RAGE"]);
        QI_LOG(verbose, "Building look-up spline");
        int num_entries = 100;
        Eigen::ArrayXd T1_values = Eigen::ArrayXd::LinSpaced(num_entries, 0.25, 4.0);
        Eigen::ArrayXd MP2_values(num_entries);
        for (int i = 0; i < num_entries; i++) {
            const auto sig = One_MP2RAGE(1., T1_values[i], 1., mp2rage_sequence);
            const float mp2 = MP2Contrast(sig[0], sig[1]);
            QI_DB(mp2)
            if ((i > 0) && (mp2 > MP2_values[i - 1])) {
                num_entries = i;
                break;
            } else {
                MP2_values[i] = mp2;
            }
        }
        QI_LOG(verbose, "There are " << num_entries << " entries in lookup table");
        QI::SplineInterpolator mp2_to_t1(MP2_values.head(num_entries), T1_values.head(num_entries));
        if (beta) {
            QI_LOG(verbose, "Recalculating unregularised MP2 image");
            MP2Filter->SetFunctor( [&](const std::complex<float> &p1, const std::complex<float> &p2) { return MP2Contrast(p1, p2, 0.0); });
            MP2Filter->Update();
        }
        using UnaryFilter = itk::UnaryGeneratorImageFilter<QI::VolumeF, QI::VolumeF>;
        auto T1LookupFilter = UnaryFilter::New();
        T1LookupFilter->SetInput(MP2Filter->GetOutput());
        auto lookup = [&](const float &p) { return mp2_to_t1(p); };
        T1LookupFilter->SetFunctor(lookup);
        QI_LOG(verbose, "Calculating T1");
        T1LookupFilter->Update();
        QI::WriteImage(T1LookupFilter->GetOutput(), out_prefix + "_MP2_T1" + QI::OutExt(), verbose);
    }
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
