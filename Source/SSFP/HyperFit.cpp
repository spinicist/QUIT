/*
 *  HyperFit.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "HyperFit.h"
#include <Eigen/Eigenvalues>
#include "Fit.h"

namespace QI {

typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

void SemiaxesToHoff(const double A, const double B, const double c,
                    double &G, double &a, double &b) {
    b = (-c*A + sqrt(c*c*A*A - (c*c + B*B)*(A*A - B*B)))/(c*c + B*B);
    a = B / (b*B + c*sqrt(1-b*b));
    G = c*(1 - b*b)/(1 - a*b);
}

Eigen::MatrixXd HyperS(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) {
    Eigen::Matrix<double, Eigen::Dynamic, 6> D(x.rows(), 6);
    D.col(0) = x*x;
    D.col(1) = 2*x*y;
    D.col(2) = y*y;
    D.col(3) = 2*x;
    D.col(4) = 2*y;
    D.col(5).setConstant(1);
    return D.transpose() * D;
}

Eigen::MatrixXd FitzC() {
    Matrix6d C = Matrix6d::Zero();
    // Fitgibbon et al
    C(0,2) = -2; C(1,1) = 1; C(2,0) = -2;
    return C;
}

Eigen::MatrixXd HyperC(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) {
    Matrix6d C = Matrix6d::Zero();
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

QI::FitReturnType HyperFit::fit(const std::vector<Eigen::ArrayXcd> &inputs,
                                const Eigen::ArrayXd &/* Unused */, QI_ARRAYN(OutputType, EllipseModel::NV) &outputs,
                                ResidualType &/*Unused*/, std::vector<Eigen::ArrayXcd> &/*Unused*/, FlagType &/*Unused*/,
                                const int /* Unused */) const
{
    const double scale = inputs[0].abs().maxCoeff();
    const Eigen::ArrayXcd data = inputs[0] / scale;
    Eigen::ArrayXd x = inputs[0].real();
    Eigen::ArrayXd y = inputs[0].imag();
    
    Eigen::MatrixXd S = HyperS(x, y);
    Matrix6d C = HyperC(x, y);
    
    // Note S and C are swapped so we can use GES
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix6d> solver(C, S);
    Vector6d Z;
    if (fabs(solver.eigenvalues()[5]) > fabs(solver.eigenvalues()[0]))
        Z = solver.eigenvectors().col(5);
    else
        Z = solver.eigenvectors().col(0);

    const double dsc=(Z[1]*Z[1]-Z[0]*Z[2]);
    const double xc = (Z[2]*Z[3]-Z[1]*Z[4])/dsc;
    const double yc = (Z[0]*Z[4]-Z[1]*Z[3])/dsc;
    const double theta_te = atan2(yc,xc);
    const double num = 2*(Z[0]*(Z[4]*Z[4])+Z[2]*(Z[3]*Z[3])+Z[5]*(Z[1]*Z[1])-2*Z[1]*Z[3]*Z[4]-Z[0]*Z[2]*Z[5]);
    double A = sqrt(num/(dsc*(sqrt((Z[0]-Z[2])*(Z[0]-Z[2]) + 4*Z[1]*Z[1])-(Z[0]+Z[2]))));
    double B = sqrt(num/(dsc*(-sqrt((Z[0]-Z[2])*(Z[0]-Z[2]) + 4*Z[1]*Z[1])-(Z[0]+Z[2]))));
    if (A > B) {
        std::swap(A, B);
    }
    double G, a, b;
    double c = sqrt(xc*xc+yc*yc);
    SemiaxesToHoff(A, B, c, G, a, b);

    /* Calculate theta_tr, i.e. drop RF phase, eddy currents etc. */
    /* First, center, rotate back to vertical and get 't' parameter */
    const Eigen::ArrayXcd vert = data / std::polar(1.0, theta_te);
    const Eigen::ArrayXd ct = (vert.real() - c) / A;
    const Eigen::VectorXd rhs = (ct - b) / (b*ct - 1);
    Eigen::MatrixXd lhs(rhs.rows(), 2);
    lhs.col(0) = cos(sequence->PhaseInc);
    lhs.col(1) = sin(sequence->PhaseInc);
    const Eigen::VectorXd K = QI::LeastSquares(lhs, rhs);
    const double theta_0 = atan2(K[1], K[0]);
    const double psi_0 = std::arg(std::polar(1.0, theta_te) / std::polar(1.0, theta_0/2));
    outputs << G * scale, a, b, theta_0  / (2*M_PI*sequence->TR), psi_0;
    return std::make_tuple(true, "");
}

} // End namespace QI