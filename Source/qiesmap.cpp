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

#include <getopt.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "QI/Util.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

class ESAlgo : public Algorithm<complex<double>> {
protected:
    size_t m_size = 0;
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 6; }
    const vector<string> & names() const {
        static vector<string> _names = {"M", "T1", "T2", "th", "a", "b"};
        return _names;
    }
    size_t dataSize() const override { return m_size; }
    void setSize(const size_t s) { m_size = s; }
    virtual TArray defaultConsts() override {
        // B1
        TArray def = TArray::Ones(0);
        return def;
    }
    void apply(const TInput &data, const TArray &inputs, TArray &outputs, TArray &resids, TIterations &its) const override
    {
        typedef Matrix<double, 6, 6> Matrix6d;
        typedef Matrix<double, 6, 1> Vector6d;
        double scale = data.abs().maxCoeff();
        ArrayXd x = data.real() / scale;
        ArrayXd y = data.imag() / scale;
        
        Matrix<double, Dynamic, 6> D(data.rows(), 6);
        D.col(0) = x*x;
        D.col(1) = x*y;
        D.col(2) = y*y;
        D.col(3) = x;
        D.col(4) = y;
        D.col(5).setConstant(1);
        Matrix6d S = D.transpose() * D;
        Matrix6d C = Matrix6d::Zero();
        C(0,2) = -2; C(1,1) = 1; C(2,0) = -2;
        
        typedef GeneralizedEigenSolver<Matrix6d> ESolver;
        ESolver es(S, C);
        ArrayXd eVals = es.eigenvalues().real();
        ArrayXd Z;
        for (int i = 0; i < 6; i++) {
            if (isfinite(eVals(i)) && (eVals(i) < 0)) {
                Z = es.eigenvectors().col(i);
                break;
            }
        }
        
        const double za = Z[0];
        const double zb = Z[1]/2;
        const double zc = Z[2];
        const double zd = Z[3]/2;
        const double zf = Z[4]/2;
        const double zg = Z[5];
        const double dsc=(zb*zb-za*zc);
        const double xc = (zc*zd-zb*zf)/dsc;
        const double yc = (za*zf-zb*zd)/dsc;
        const double th = atan2(yc,xc);
        double A = sqrt((2*(za*(zf*zf)+zc*(zd*zd)+zg*(zb*zb)-2*zb*zd*zf-za*zc*zg))/(dsc*(sqrt((za-zc)*(za-zc) + 4*zb*zb)-(za+zc))));
        double B = sqrt((2*(za*(zf*zf)+zc*(zd*zd)+zg*(zb*zb)-2*zb*zd*zf-za*zc*zg))/(dsc*(-sqrt((za-zc)*(za-zc) + 4*zb*zb)-(za+zc))));
        if (A > B) {
            std::swap(A, B);
        }
        cout << "Z " << Z.transpose() << endl;
        cout << "dsc " << dsc << " xc " << xc << " yc " << yc << " A " << A << " B " << B << endl;
        const double c = sqrt(xc*xc+yc*xc);
        const double b = (-2*c*A + sqrt(pow(2*c*A,2)-4*(c*c+B*B)*(A*A-B*B)))/(2*(c*c+B*B));
        const double a = B / (b*B + c*sqrt(1-b*b));
        const double M = scale*c*(1-b*b)/(1-a*b);
        const double TR = 0.0065;
        const double FA = (25*M_PI/180.);
        const double T1 = -TR / log((a*(1+cos(FA)-a*b*cos(FA))-b)/(a*(1+cos(FA)-a*b)-b*cos(FA)));
        const double T2 = -TR / log(a);
        
        outputs[0] = M;
        outputs[1] = T1;
        outputs[2] = T2;
        outputs[3] = th;
        outputs[4] = a;
        outputs[5] = b;
    }
};

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qiesmap [options] input \n\
\n\
A utility for calculating T1,T2,PD and f0 maps from SSFP data.\n\
Input must be a single complex image with at least 6 phase increments.\n\
\n\
Options:\n\
    --help, -h       : Print this message.\n\
    --verbose, -v    : Print more information.\n\
    --out, -o path   : Specify an output prefix.\n\
    --mask, -m file  : Mask input with specified file.\n\
    --threads, -T N  : Use N threads (default=hardware limit).\n"
};

bool verbose = false;
static size_t num_threads = 4;
static string outPrefix;
const struct option long_opts[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"out", required_argument, 0, 'o'},
    {"mask", required_argument, 0, 'm'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
const char *short_opts = "hvo:m:T:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::VolumeF::Pointer mask = ITK_NULLPTR;
    QI::VolumeF::Pointer B1   = ITK_NULLPTR;

    shared_ptr<ESAlgo> algo = make_shared<ESAlgo>();
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'm':
                if (verbose) cout << "Opening mask file " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
            case 'o':
                outPrefix = optarg;
                if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
                break;
            case 'b':
                if (verbose) cout << "Opening B1 file: " << optarg << endl;
                B1 = QI::ReadImage(optarg);
                break;
            case 'T':
                num_threads = stoi(optarg);
                if (num_threads == 0)
                    num_threads = std::thread::hardware_concurrency();
                break;
            case 'h':
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 1) {
        cout << "Incorrect number of arguments." << endl << usage << endl;
        return EXIT_FAILURE;
    }

    string inputFilename = argv[optind++];
    if (verbose) cout << "Opening file: " << inputFilename << endl;
    auto data = QI::ReadVectorImage<complex<float>>(inputFilename);
    auto apply = itk::ApplyAlgorithmFilter<ESAlgo, complex<float>>::New();
    algo->setSize(data->GetNumberOfComponentsPerPixel());
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(num_threads);
    apply->SetInput(0, data);
    if (mask)
        apply->SetMask(mask);
    /*if (B1)
        apply->SetConst(0, B1);*/
    if (verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    outPrefix = outPrefix + "ES_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        QI::WriteImage(apply->GetOutput(i), outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
