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

#include <algorithm>
#include <fstream>
#include <string>

#include <Eigen/Dense>

#include "MeanImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkTileImageFilter.h"

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"

namespace QI {
template <typename T> bool Read(const std::string &s, T &val) {
    std::istringstream stream(s);
    if (!(stream >> val)) {
        QI::Fail("Failed to parse input: {}", s);
    }
    return true;
}

template <typename T> bool Read(std::istream &in, T &val) {
    std::string line;
    // Ignore comment lines. Use shell script convention
    while (in.peek() == '#') {
        if (!std::getline(in, line))
            QI::Fail("Failed to read input.");
    }
    if (!std::getline(in, line)) {
        QI::Fail("Failed to read input. Last line was: {}", line);
    }
    return Read(line, val);
}

template <typename Scalar>
void ReadArray(const std::string &s, Eigen::Array<Scalar, Eigen::Dynamic, 1> &array) {
    std::istringstream  stream(s);
    std::vector<Scalar> vals;

    Scalar temp;
    while (stream >> temp) {
        vals.push_back(temp);
    }

    array = Eigen::Array<Scalar, Eigen::Dynamic, 1>(vals.size());
    for (size_t i = 0; i < vals.size(); i++) {
        array[i] = vals[i];
    }
}

template <typename Scalar>
void ReadArray(std::istream &in, Eigen::Array<Scalar, Eigen::Dynamic, 1> &array) {
    std::string line;
    // Ignore comment lines. Use shell script convention
    while (in.peek() == '#') {
        if (!std::getline(in, line))
            QI::Fail("Failed to read line: {}", line);
    }
    if (!std::getline(in, line)) {
        QI::Fail("Failed to read line: {}", line);
    }
    ReadArray(line, array);
}

template <typename Scalar>
void ReadArray(const std::string &s, Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &array) {
    std::istringstream  stream(s);
    std::vector<Scalar> vals;

    Scalar temp;
    while (stream >> temp) {
        vals.push_back(temp);
    }

    array = Eigen::Array<Scalar, Eigen::Dynamic, 1>(vals.size());
    for (int i = 0; i < vals.size(); i++) {
        array[i] = vals[i];
    }
}

template <typename Scalar>
void ReadArray(std::istream &in, Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &array) {
    std::vector<Eigen::Array<Scalar, Eigen::Dynamic, 1>> rows;
    std::string                                          line;

    Eigen::Array<Scalar, Eigen::Dynamic, 1> row;
    while (std::getline(in, line) && line != "") {
        ReadArray(line, row);
        rows.push_back(row);
    }

    array = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>(rows.size(), rows.at(0).size());
    for (size_t i = 0; i < rows.size(); i++) {
        if (rows.at(i).size() == array.cols()) {
            array.row(i) = rows.at(i);
        } else {
            QI::Fail("Inconsistent row sizes in input");
        }
    }
}

// Utility function to read a file into an array, ignoring comment lines
Eigen::ArrayXXd ReadArrayFile(const std::string &path) {
    std::ifstream matrix_file(path);
    if (!matrix_file) {
        QI::Fail("Failed to open matrix file: {}", path);
    }
    while (matrix_file && (matrix_file.peek() == '/')) {
        matrix_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    Eigen::ArrayXXd matrix;
    QI::ReadArray(matrix_file, matrix);
    return matrix;
}

} // End namespace QI

class ContrastsFilter : public itk::ImageToImageFilter<QI::VectorVolumeF, QI::VolumeF> {
  public:
    /** Standard class typedefs. */
    typedef ContrastsFilter                                    Self;
    typedef ImageToImageFilter<QI::VectorVolumeF, QI::VolumeF> Superclass;
    typedef itk::SmartPointer<Self>                            Pointer;
    typedef typename QI::VolumeF::RegionType                   RegionType;

    itkNewMacro(Self)
    itkTypeMacro(Self, Superclass);

    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto input     = this->GetInput(0);
        auto region    = input->GetLargestPossibleRegion();
        auto spacing   = input->GetSpacing();
        auto origin    = input->GetOrigin();
        auto direction = input->GetDirection();
        for (int i = 0; i < m_mat.rows(); i++) {
            auto op = this->GetOutput(i);
            op->SetRegions(region);
            op->SetSpacing(spacing);
            op->SetOrigin(origin);
            op->SetDirection(direction);
            op->Allocate(true);
        }
    }

    void SetMatrix(const Eigen::MatrixXd &d, const Eigen::MatrixXd &c, const bool s = false) {
        m_scale = s;
        m_mat   = c * ((d.transpose() * d).inverse()) * d.transpose();
        this->SetNumberOfRequiredOutputs(m_mat.rows());
        for (int i = 0; i < m_mat.rows(); i++) {
            this->SetNthOutput(i, this->MakeOutput(i));
        }
    }

  protected:
    bool            m_scale = false;
    Eigen::MatrixXd m_mat;

    ContrastsFilter() { this->SetNumberOfRequiredInputs(1); }
    ~ContrastsFilter() {}

    void DynamicThreadedGenerateData(const RegionType &region) ITK_OVERRIDE {
        const auto                                       input_image = this->GetInput(0);
        itk::ImageRegionConstIterator<QI::VectorVolumeF> input_iter(input_image, region);
        input_iter.GoToBegin();
        std::vector<itk::ImageRegionIterator<QI::VolumeF>> out_iters(m_mat.rows());
        for (int i = 0; i < m_mat.rows(); i++) {
            out_iters.at(i) = itk::ImageRegionIterator<QI::VolumeF>(this->GetOutput(i), region);
            out_iters.at(i).GoToBegin();
        }

        while (!input_iter.IsAtEnd()) {
            const auto input_vec = input_iter.Get();

            Eigen::Map<const Eigen::VectorXf> indata(input_vec.GetDataPointer(), input_vec.Size());
            Eigen::VectorXd                   c = m_mat * indata.cast<double>();
            if (m_scale) {
                c /= indata.mean();
            }
            for (int i = 0; i < m_mat.rows(); i++) {
                out_iters[i].Set(c[i]);
                ++out_iters[i];
            }
            ++input_iter;
        }
    }

  private:
    ContrastsFilter(const Self &); // purposely not implemented
    void operator=(const Self &);  // purposely not implemented
};

/*
 * Main
 */
int glm_contrasts_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(
        parser, "IMAGE", "The combined image file from qi_glmsetup");
    args::Positional<std::string> design_path(
        parser, "DESIGN", "GLM Design matrix from qi_glmsetup");
    args::Positional<std::string> contrasts_path(
        parser, "CONTRASTS", "Contrasts matrix from qi_glmsetup");
    args::Flag fraction(
        parser, "FRACTION", "Output contrasts as fraction of grand mean", {'F', "frac"});
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    parser.Parse();

    QI::Log(verbose, "Reading input file {}", QI::CheckPos(input_path));
    QI::VectorVolumeF::Pointer merged =
        QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    QI::Log(verbose, "Reading design matrix {}", QI::CheckPos(design_path));
    Eigen::ArrayXXd design_matrix = QI::ReadArrayFile(QI::CheckPos(design_path));
    QI::Log(verbose, "Reading contrasts file {}", QI::CheckPos(contrasts_path));
    Eigen::ArrayXXd contrasts = QI::ReadArrayFile(QI::CheckPos(contrasts_path));
    if (design_matrix.rows() != merged->GetNumberOfComponentsPerPixel()) {
        QI::Fail(
            "Number of rows in design matrix ({}) does not match number of volumes in image ({})",
            design_matrix.rows(),
            merged->GetNumberOfComponentsPerPixel());
    }
    if (design_matrix.cols() != contrasts.cols()) {
        QI::Fail("Number of columns in design matrix ({}) does not match contrasts ({})",
                 design_matrix.cols(),
                 contrasts.cols());
    }

    auto con_filter = ContrastsFilter::New();
    con_filter->SetMatrix(design_matrix, contrasts, fraction);
    con_filter->SetInput(0, merged);
    QI::Log(verbose, "Calculating contrasts");
    con_filter->Update();
    for (int c = 0; c < contrasts.rows(); c++) {
        QI::Log(verbose, "Writing contrast {}", c + 1);
        QI::WriteImage(con_filter->GetOutput(c),
                       outarg.Get() + "con" + std::to_string(c + 1) + QI::OutExt(),
                       verbose);
    }
    return EXIT_SUCCESS;
}
