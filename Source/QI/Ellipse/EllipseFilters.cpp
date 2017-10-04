/*
 *  Ellipse.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Ellipse/EllipseFilters.h"
/*#include "QI/Ellipse/EllipseHelpers.h"
#include "QI/Banding.h"
#include "QI/GoldenSection.h"*/

namespace QI {

/*void ConstrainedEllipse::fit(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y, const Eigen::Vector2d p, const Eigen::Vector2d q, 
             Eigen::Vector2d gammaBound, Eigen::Matrix2d &A, Eigen::Vector2d &x_c, double &g) const {
    const double epsFmin = 1e-8, epsRootsImag = 1e-6, epsRootsMin = 1e-3; // Threshold Values
    Eigen::Matrix3d B; B << 0,  0, 2,
                            0, -1, 0,
                            2,  0, 0;                                     // Constraint Matrix
    Eigen::PartialPivLU<Eigen::Matrix3d> BLU(B);
    const size_t n = x.rows();
    Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(n,n) - Eigen::MatrixXd::Ones(n,n)/n; // Centering Matrix
    Eigen::MatrixXd D0(n, 3); D0.col(0) = (x-p[0]).square(); D0.col(1) = (x-p[0])*(y-p[1]); D0.col(2) = (y-p[1]).square();
    Eigen::MatrixXd D1(n, 3); D1.col(0) = 2*q[0]*(x-p[0]); D1.col(1) = q[0]*(y-p[1])+q[1]*(x-p[0]); D1.col(2) = 2*q[1]*(y-p[1]);
    Eigen::Matrix3d C0 = D0.transpose()*Z*D0;
    Eigen::Matrix3d C1 = -D0.transpose()*Z*D1 - D1.transpose()*Z*D0;
    Eigen::Matrix3d C2 = D1.transpose()*Z*D1;
    Eigen::Matrix3d BC0 = BLU.solve(C0);
    Eigen::Matrix3d BC1 = BLU.solve(C1);
    Eigen::Matrix3d BC2 = BLU.solve(C2);

    std::function<double(double)> getLambdaOfGamma = [&] (double gamma) {
        Eigen::Matrix3d M = BC0 + gamma*BC1 + (gamma*gamma)*BC2;
        Eigen::EigenSolver<Eigen::Matrix3d> eig(M);
        double lambda = fabs(eig.eigenvalues().real().maxCoeff());
        return lambda;
    };

    bool isGlobalMinimumFound = false;
    int its = 0;
    while (!isGlobalMinimumFound) {
        // GetLamdaOfGamma
        double gamma = GoldenSectionSearch(getLambdaOfGamma, gammaBound[0], gammaBound[1], epsFmin);
        double lambda = getLambdaOfGamma(gamma);
        //cout << "gamma " << gamma << " lambda " << lambda << endl;
        Eigen::Matrix3d P; P << 1, 0, q[1]*q[1],
                                0, 1, -2*q[0]*q[1],
                                0, 0, q[0]*q[0];
        Eigen::Matrix3d C2t = P.transpose()*C2*P;
        Eigen::Matrix3d C1t = P.transpose()*C1*P;
        Eigen::Matrix3d C0t = P.transpose()*C0*P;
        Eigen::Matrix3d Bt = P.transpose()*B*P;
        C0t = C0t - lambda*Bt;
        Eigen::Matrix2d C2b = C2t.topLeftCorner<2,2>()-C1t.block<2,1>(0,2)*C1t.block<1,2>(2,0)/C0t.coeff(2,2);
        Eigen::Matrix2d C1b = C1t.topLeftCorner<2,2>()-C1t.block<2,1>(0,2)*C0t.block<1,2>(2,0)/C0t.coeff(2,2) - C0t.block<2,1>(0,2)*C1t.block<1,2>(2,0)/C0t.coeff(2,2);
        Eigen::Matrix2d C0b = C0t.topLeftCorner<2,2>()-C0t.block<2,1>(0,2)*C0t.block<1,2>(2,0)/C0t.coeff(2,2);
        Eigen::Matrix4d rootsA; rootsA << Eigen::Matrix2d::Zero(), Eigen::Matrix2d::Identity(),
                                            C0b, C1b;
        Eigen::Matrix4d rootsB; rootsB << Eigen::Matrix2d::Identity(), Eigen::Matrix2d::Zero(),
                                            Eigen::Matrix2d::Zero(), -C2b;
        Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> roots(rootsA, rootsB);
        Eigen::Vector4cd rootsVec = roots.eigenvalues();
        //cout << "rootsVec " << rootsVec.transpose() << endl;
        auto realRootsIdx = (rootsVec.imag().array().abs() < epsRootsImag) && ((rootsVec.array() - gamma).abs() > epsRootsMin);
        //cout << "realRootsIdx " << realRootsIdx.transpose() << endl;
        std::vector<double> realRoots; for (int i = 0; i < 4; i++) { if (realRootsIdx[i]) realRoots.push_back(real(rootsVec[i])); }
        std::sort(realRoots.begin(), realRoots.end());
        // cout << realRoots.size() << " : "; for (int i = 0; i < realRoots.size(); i++) { cout << realRoots[i] << " "; } cout << endl;
        if (realRoots.size() == 0) {
            isGlobalMinimumFound = true;
        } else if (realRoots.size() == 1) {
            // 1 real root. Global minimum is outside the given search range
            if (realRoots[0] < gammaBound[0]) {
                gammaBound[1] = gammaBound[0]; gammaBound[0] = realRoots[0];
            } else {
                gammaBound[0] = gammaBound[1]; gammaBound[1] = realRoots[0];
            }
        } else if (realRoots.size() == 2) {
            if (((realRoots[0] < gammaBound[0]) && (realRoots[1] > gammaBound[1])) || // 2 roots outside bounds
                ((realRoots[0] > gammaBound[0]) && (realRoots[1] < gammaBound[1]))) { // 2 roots inside bounds
                gammaBound[0] = realRoots[0]; gammaBound[1] = realRoots[1];
            } else if (realRoots[0] < gammaBound[0]) { // 1 root below bounds
                gammaBound[1] = realRoots[1];
            } else {
                gammaBound[0] = realRoots[0];
            }
            if (gamma > gammaBound[0] && gamma < gammaBound[1]) {
                // This shouldn't happen, but might due to noise. In which case claim it is the global minimum
                isGlobalMinimumFound = true;
            }
        } else { // 3 roots
            bool allOutside = true;
            std::vector<float> newBounds;
            for (int i = 0; i < 3; i++) {
                if (realRoots[i] > gammaBound[0] && realRoots[i] < gammaBound[1]) {
                    allOutside = false;
                    newBounds.push_back(realRoots[i]);
                }
            }
            if (allOutside) {
                isGlobalMinimumFound = true;
            } else {
                gammaBound[0] = std::min(newBounds[0], newBounds[1]);
                gammaBound[1] = std::max(newBounds[0], newBounds[1]);
            }
        }
        if (gammaBound[0] == gammaBound[1])
            isGlobalMinimumFound = true;
        //cout << "new gammaBound " << gammaBound.transpose() << endl;
        Eigen::Matrix3d M = BC0 + gamma*BC1 + gamma*gamma*BC2;
        Eigen::EigenSolver<Eigen::Matrix3d> eig(M);
        int eigIdx; eig.eigenvalues().real().maxCoeff(&eigIdx);
        //cout << "M eigenvalues " << eig.eigenvalues().transpose() << " idx " << eigIdx << endl;
        Eigen::Vector3d theta = eig.eigenvectors().col(eigIdx).real();
        // First part is a hacked sign function
        //cout << "theta " << theta.transpose() << endl;
        theta = ((theta[0] > 0) - (theta[0] < 0))*theta/sqrt(theta.transpose()*B*theta); // Normalization
        //cout << "theta n" << theta.transpose() << endl;
        double h_opt = (-1.0/n * Eigen::VectorXd::Ones(n).transpose()) * (D0 - gamma*D1)*theta;
        A << theta[0], theta[1]/2,
                theta[1]/2, theta[2];
        x_c = p + gamma*q;
        g = h_opt - gamma*gamma*q.transpose()*A*q;
        its++;
        if (its > 10)
            break;
    }
}*/

} // End namespace QI