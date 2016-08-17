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
#include <array>
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
    T2 = -TR / log(a);
    T1 = -TR / (log(a-b + a*cosf*(1.-a*b)) - log(a*(1.-a*b) + (a-b)*cosf));
    //cout << "TR " << TR << " a " << a << " cla " << clampa << " b " << b << " T1 " << T1 << endl;
    M = (scale/sqrt(a))*(1-b*b)/(1-a*b);
    df0 = th / (M_PI*TR); // Missing factor of 2 due to TE not TR
}

class ESAlgo : public QI::ApplyVectorXFVectorF::Algorithm {
protected:
    bool m_reorderPhase = false;
    shared_ptr<QI::SSFP_GS> m_sequence = nullptr;
    size_t m_pincs = 6;
    TOutput m_zero;
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 6; }
    size_t dataSize() const override { return m_sequence->size() * m_pincs; }
    size_t outputSize(const int i) const override { return m_sequence->flip().rows(); }
    void setReorderPhase(const bool p) { m_reorderPhase = p; }
    void SetSequence(const shared_ptr<QI::SSFP_GS> &s, const size_t pincs) {
        m_sequence = s;
        m_pincs = pincs;
        m_zero = TOutput(m_sequence->flip().rows());
        m_zero.Fill(0.);
    }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 1.0f); // B1
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override { return m_zero; }
    const vector<string> & names() const {
        static vector<string> _names = {"M", "T1", "T2", "f0", "a", "b"};
        return _names;
    }
    
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

    std::array<float, 6> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &indata, const double B1) const {
        typedef Matrix<double, 6, 6> Matrix6d;
        typedef Matrix<double, 6, 1> Vector6d;
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
        std::array<float, 6> outputs;
        EllipseToMRI(a, b, c*scale, th, m_sequence->TR(), B1 * m_sequence->flip()[0],
                     outputs[0], outputs[1], outputs[2], outputs[3]);
        outputs[4] = a;
        outputs[5] = b;
        return outputs;
    }

    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        const double B1 = consts[0];
        size_t phase_stride = 1;
        size_t flip_stride = m_sequence->flip().rows();
        if (m_reorderPhase)
            std::swap(phase_stride, flip_stride);
        for (int f = 0; f < m_sequence->flip().rows(); f++) {
            const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> vf(inputs[0].GetDataPointer() + f*flip_stride, m_pincs, Eigen::InnerStride<>(phase_stride));
            std::array<float, 6> tempOutputs = this->applyFlip(vf, B1);
            for (int o = 0; o < 6; o++) {
                outputs[o][f] = tempOutputs[o];
            }
        }
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts("Usage is: qiesmap [options] input_file\n\nA utility for calculating T1,T2,PD and f0 maps from SSFP data.\nInput must be a single complex image with at least 6 phase increments.");
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::Option<int> ph_incs(6,'\0',"ph_incs","Number of phase increments (default is 6).", opts);
    QI::Switch ph_order('\0',"ph_order","Data order is phase, then flip-angle (default opposite).", opts);
    QI::EnumOption algorithm("lwnb",'l','a',"algo","Choose algorithm (f/h/c)", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Prefix output filenames", opts);
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
    auto apply = QI::ApplyVectorXFVectorF::New();
    algo->SetSequence(seq, *ph_incs);
    algo->setReorderPhase(*ph_order);
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
        QI::WriteVectorImage(apply->GetOutput(i), *outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
