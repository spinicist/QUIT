/*
 *  qiesmap.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "QI/Util.h"
#include "QI/Option.h"
#include "QI/Sequences/SteadyStateSequence.cpp"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

// Helper Functions
void SemiaxesToHoff(const double A, const double B, const double c,
                    double &a, double &b) {
    b = (-c*A + sqrt(c*c*A*A - (c*c + B*B)*(A*A - B*B)))/(c*c + B*B);
    a = B / (b*B + c*sqrt(1-b*b));
}


void EllipseToMRI(const double a, const double b, const double scale, const double th, const double TR, const double flip,
                  float &M, float &T1, float &T2, float &df0) {
    const double cosf = cos(flip);
    const double upper = exp(-TR / 4.3);
    const double clampa = QI::clamp(a, 0.0, upper);
    T2 = -TR / log(clampa);
    T1 = -TR / (log(clampa-b + (1.-clampa*b)*clampa*cosf) - log(clampa*(1.-clampa*b) + (clampa-b)*cosf));
    
    M = (scale/sqrt(clampa))*(1-b*b)/(1-clampa*b);
    df0 = th / (2.*M_PI*TR);
}

class ESAlgo : public QI::ApplyXF::Algorithm {
protected:
    size_t m_size = 0;
    shared_ptr<QI::SSFP_GS> m_sequence = nullptr;
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 6; }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }
    const vector<string> & names() const {
        static vector<string> _names = {"M", "T1", "T2", "th", "a", "b"};
        return _names;
    }
    size_t dataSize() const override { return m_size; }
    void setSize(const size_t s) { m_size = s; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 1.0f); // B1
        return def;
    }
    void SetSequence(const shared_ptr<QI::SSFP_GS> &s) { m_sequence = s;}
    
    MatrixXd buildS(const ArrayXd &x, const ArrayXd &y) const {
        Matrix<double, Dynamic, 6> D(x.rows(), 6);
        D.col(0) = x*x;
        D.col(1) = 2*x*y;
        D.col(2) = y*y;
        D.col(3) = 2*x;
        D.col(4) = 2*y;
        D.col(5).setConstant(1);
        return D.transpose() * D;
    }
    
    MatrixXd fitzC() const {
        typedef Matrix<double, 6, 6> Matrix6d;
        Matrix6d C = Matrix6d::Zero();
        // FitZ[5]ibbon et al
        C(0,2) = -2; C(1,1) = 1; C(2,0) = -2;
        return C;
    }
    
    MatrixXd hyperC(const ArrayXd &x, const ArrayXd &y) const {
        typedef Matrix<double, 6, 6> Matrix6d;
        Matrix6d C = Matrix6d::Zero();
        // FitZ[5]ibbon et al
        //C(0,2) = -2; C(1,1) = 1; C(2,0) = -2;
        
        // Hyper Ellipse
        const double N = x.cols();
        const double xc = x.sum() / N;
        const double yc = y.sum() / N;
        const double sx = x.square().sum() / N;
        const double sy = y.square().sum() / N;
        const double xy = (x * y).sum() / N; 
        
        C << 6*sx, 6*xy, sx+sy, 6*xc, 2*yc, 1,
             6*xy, 4*(sx+sy), 6*xy, 4*yc, 4*xc, 0,
             sx + sy, 6*xy, 6*sy, 2*xc, 6*yc, 1,
             6*xc, 4*yc, 2*xc, 4, 0, 0,
             2*yc, 4*xc, 6*yc, 0, 4, 0,
             1, 0, 1, 0, 0, 0;
        
        return C;
    }
    
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        typedef Matrix<double, 6, 6> Matrix6d;
        typedef Matrix<double, 6, 1> Vector6d;
        const double B1 = consts[0];
        Eigen::Map<const Eigen::ArrayXcf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        ArrayXcd data = indata.cast<complex<double>>();
        const double scale = data.abs().maxCoeff();
        ArrayXd x = data.real() / scale;
        ArrayXd y = data.imag() / scale;
        
        MatrixXd S = buildS(x, y);
        Matrix6d C = hyperC(x, y);
        
        // Note A and B are swapped so we can use GES
        GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(C, S);
        ArrayXd Z;
        if (fabs(solver.eigenvalues()[5]) > fabs(solver.eigenvalues()[0]))
            Z = solver.eigenvectors().col(5);
        else
            Z = solver.eigenvectors().col(0);

        const double dsc=(Z[1]*Z[1]-Z[0]*Z[2]);
        const double xc = (Z[2]*Z[3]-Z[1]*Z[4])/dsc;
        const double yc = (Z[0]*Z[4]-Z[1]*Z[3])/dsc;
        const double th = atan2(yc,xc);
        const double num = 2*(Z[0]*(Z[4]*Z[4])+Z[2]*(Z[3]*Z[3])+Z[5]*(Z[1]*Z[1])-2*Z[1]*Z[3]*Z[4]-Z[0]*Z[2]*Z[5]);
        double A = sqrt(num/(dsc*(sqrt((Z[0]-Z[2])*(Z[0]-Z[2]) + 4*Z[1]*Z[1])-(Z[0]+Z[2]))));
        double B = sqrt(num/(dsc*(-sqrt((Z[0]-Z[2])*(Z[0]-Z[2]) + 4*Z[1]*Z[1])-(Z[0]+Z[2]))));
        if (A > B) {
            std::swap(A, B);
        }
        double a, b;
        double c = sqrt(xc*xc+yc*yc);
        SemiaxesToHoff(A, B, c, a, b);
        EllipseToMRI(a, b, c*scale, th, m_sequence->TR(), B1 * m_sequence->flip()[0],
                     outputs[0], outputs[1], outputs[2], outputs[3]);
        outputs[4] = a;
        outputs[5] = b;
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts("Usage is: qiesmap [options] input_file\n\nA utility for calculating T1,T2,PD and f0 maps from SSFP data.\nInput must be a single complex image with at least 6 phase increments.");
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Prefix output filenames", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::EnumOption algorithm("lwnb",'l','a',"algo","Choose algorithm (f/h/c)", opts);
    QI::Switch suppress('n',"no-prompt","Suppress input prompts", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 1) {
        std::cerr << opts << std::endl;
        std::cerr << "No input filename specified." << std::endl;
        return EXIT_FAILURE;
    }
    shared_ptr<ESAlgo> algo = make_shared<ESAlgo>();
    if (*verbose) cout << "Opening file: " << nonopts[0] << endl;
    auto data = QI::ReadVectorImage<complex<float>>(nonopts[0]);
    shared_ptr<QI::SSFP_GS> seq = make_shared<QI::SSFP_GS>(cin, !*suppress);
    if (*verbose) cout << *seq;
    auto apply = QI::ApplyXF::New();
    algo->setSize(data->GetNumberOfComponentsPerPixel());
    algo->SetSequence(seq);
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(*num_threads);
    apply->SetInput(0, data);
    apply->SetMask(*mask);
    apply->SetConst(0, *B1);
    if (*verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (*verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    *outPrefix = *outPrefix + "ES_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        QI::WriteImage(apply->GetOutput(i), *outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
